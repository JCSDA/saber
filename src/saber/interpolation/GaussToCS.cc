/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include "atlas/field.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

#include "saber/interpolation/GaussToCS.h"
#include "saber/interpolation/Rescaling.h"

using atlas::grid::detail::partitioner::TransPartitioner;
using atlas::grid::detail::partitioner::MatchingMeshPartitionerCubedSphere;

namespace saber {
namespace interpolation {

namespace {

int getHalo(const std::string & interpType) {
  oops::Log::trace()  << "::GaussToCS::getHalo" << std::endl;

  // Assumes that only interpolation types with "cubic" in their name require a halo of 2.
  const std::string pattern = "cubic";
  if (interpType.find(pattern) != interpType.npos) {
    return 2;
  }
  return 1;
}

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid,
                             const int halo) {
  oops::Log::trace() << "::GaussToCS::createGaussFunctionSpace" << std::endl;
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(halo));
}

// -----------------------------------------------------------------------------

auto createPointCloud(const atlas::Grid& grid,
                      const atlas::grid::Partitioner& partitioner) {
    oops::Log::trace()  << "::GaussToCS::createPointCloud" << std::endl;
    const auto distribution = atlas::grid::Distribution(grid, partitioner);

    auto lonLats = std::vector<atlas::PointXY>{};
    auto idx = atlas::gidx_t{0};
    for (const auto& lonLat : grid.lonlat()) {
      if (distribution.partition(idx++) == atlas::mpi::rank()) {
        lonLats.emplace_back(lonLat.data());
      }
    }
    return std::make_unique<atlas::functionspace::PointCloud>(
                atlas::functionspace::PointCloud(lonLats));
}

// -----------------------------------------------------------------------------

auto createInverseInterpolation(const bool initializeInverseInterpolation,
                                const bool hasSinglePE,
                                const atlas::functionspace::NodeColumns & csFunctionSpace,
                                const atlas::Grid & gaussGrid,
                                const atlas::grid::Partitioner & gaussPartitioner,
                                const std::string & interpType,
                                const bool includingVectorInterpolation) {
  oops::Log::trace()  << "::GaussToCS::createInverseInterpolation" << std::endl;
  CS2Gauss inverseInterpolation;

  if (!initializeInverseInterpolation || hasSinglePE) {
    return inverseInterpolation;
  }

  const auto config = atlas::util::Config("type", "cubedsphere");
  const atlas::grid::MatchingMeshPartitioner csMatchingPartitioner(csFunctionSpace.mesh(),
                                                                   config);
  inverseInterpolation.matchingPtcldFspace = createPointCloud(gaussGrid, csMatchingPartitioner);

  inverseInterpolation.targetPtcldFspace = createPointCloud(gaussGrid, gaussPartitioner);

  atlas::util::Config interpScalarConfig;
  interpScalarConfig.set("type", "cubedsphere-bilinear");
  interpScalarConfig.set("halo_exchange", false);

  atlas::util::Config interpConfig;
  if (includingVectorInterpolation) {
    static const auto sphericalVector =
      atlas::option::type("spherical-vector");
    interpConfig = sphericalVector | atlas::util::Config("scheme", interpScalarConfig);
  } else {
    interpConfig = interpScalarConfig;
  }

  inverseInterpolation.interpolation = atlas::Interpolation(
              interpConfig, csFunctionSpace,
              *inverseInterpolation.matchingPtcldFspace);

  inverseInterpolation.redistribution = atlas::Redistribution(
              *inverseInterpolation.matchingPtcldFspace,
              *inverseInterpolation.targetPtcldFspace);

  return inverseInterpolation;
}

// -----------------------------------------------------------------------------

/* Direct interpolation from NodeColumn cubed-sphere FunctionSpace to
 * StructuredColumns is not possible on multiple PEs.
 * Here, the interpolation is done following this route:
 * CubedSphere FunctionSpace
 * -> matching PointCloud FunctionSpace using CS partitioner
 * -> target PointCloud FunctionSpace using Gauss grid partitioner
 * -> gauss StructuredColumns FunctionSpace.
 */
void inverseInterpolateMultiplePEs(
        const oops::Variables & variables,
        const CS2Gauss & inverseInterpolation,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const bool includingVectorInterpolation,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  oops::Log::trace()  << "::GaussToCS::createInverseInterpolateMultiplePEs" << std::endl;

  srcFieldSet.haloExchange();

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = srcFieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // Interpolate from source to matching PointCloud and create srcFieldSet copy.
  atlas::FieldSet tmpSrcFieldSet;
  atlas::FieldSet matchingPtcldFset;

  for (auto & fieldname : sortedFieldNames) {
    if (variables.has(fieldname)) {
      auto matchingPtcldField =
        inverseInterpolation.matchingPtcldFspace->createField<double>(
                               atlas::option::name(fieldname) |
                               atlas::option::levels(srcFieldSet[fieldname].shape(1)));

      auto tmpSrcField =
        srcFieldSet[fieldname].functionspace().createField<double>(
                                 atlas::option::name(fieldname) |
                                 atlas::option::levels(srcFieldSet[fieldname].shape(1)));

      atlas::array::make_view<double, 2>(matchingPtcldField).assign(0.0);

      const auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[fieldname]);
      auto srcTmpView = atlas::array::make_view<double, 2>(tmpSrcField);
      for (atlas::idx_t t = 0; t < srcFieldSet[fieldname].shape(0); ++t) {
        for (atlas::idx_t k = 0; k < srcFieldSet[fieldname].shape(1); ++k) {
          srcTmpView(t, k) = srcView(t, k);
        }
      }

      matchingPtcldFset.add(matchingPtcldField);
      tmpSrcFieldSet.add(tmpSrcField);
    }
  }

  if (includingVectorInterpolation) {
    tmpSrcFieldSet["eastward_wind"].metadata().set("vector_field_name", "wind");
    tmpSrcFieldSet["northward_wind"].metadata().set("vector_field_name", "wind");
    matchingPtcldFset["eastward_wind"].metadata().set("vector_field_name", "wind");
    matchingPtcldFset["northward_wind"].metadata().set("vector_field_name", "wind");
  }

  inverseInterpolation.interpolation.execute(tmpSrcFieldSet, matchingPtcldFset);

  // Redistribute from matching PointCloud to target PointCloud
  atlas::FieldSet targetPtcldFset;
  for (auto & fieldname : sortedFieldNames) {
    if (variables.has(fieldname)) {
      auto targetPtcldField =
        inverseInterpolation.targetPtcldFspace->createField<double>(
                               atlas::option::name(fieldname) |
                               atlas::option::levels(srcFieldSet[fieldname].shape(1)));
      targetPtcldFset.add(targetPtcldField);
    }
  }

  inverseInterpolation.redistribution.execute(matchingPtcldFset, targetPtcldFset);

  // Copy from target PointCloud to gauss StructuredColumns
  for (auto & fieldname : sortedFieldNames) {
    if (variables.has(fieldname)) {
      atlas::Field gaussField =
        gaussFunctionSpace.createField<double>(
                             atlas::option::name(fieldname) |
                             atlas::option::levels(srcFieldSet[fieldname].shape(1)));
      atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(targetPtcldFset[fieldname]));
      gaussField.set_dirty();  // atlas interpolation/redistribution above produces dirty halos
      newFieldSet.add(gaussField);
    }
  }
}

// -----------------------------------------------------------------------------

void inverseInterpolateSinglePE(
        const oops::Variables & variables,
        const atlas::functionspace::NodeColumns & CSFunctionSpace,
        const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
        const atlas::Grid & gaussGrid,
        const bool includingVectorInterpolation,
        const std::string & interpType,
        const atlas::FieldSet & srcFieldSet,
        atlas::FieldSet & newFieldSet) {
  oops::Log::trace() << "::GaussToCS::createInverseInterpolateSinglePE" << std::endl;
  const auto targetMesh = atlas::MeshGenerator("structured").generate(gaussGrid);
  const auto hybridFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);

  srcFieldSet.haloExchange();

  atlas::util::Config interpScalarConfig;
  interpScalarConfig.set("type", "cubedsphere-bilinear");
  interpScalarConfig.set("adjoint", false);

  atlas::util::Config interpConfig;
  if (includingVectorInterpolation) {
    static const auto sphericalVector =
      atlas::option::type("spherical-vector");
    interpConfig = sphericalVector | atlas::util::Config("scheme", interpScalarConfig);
  } else {
    interpConfig = interpScalarConfig;
  }

  oops::Log::info() << "interpConfig inverseInterpolateSinglePE is: " << interpConfig << std::endl;
  oops::Log::info() << "CSFunctionSpace.type() is: " << CSFunctionSpace.type() << std::endl;

  const auto interp = atlas::Interpolation(interpConfig, CSFunctionSpace, hybridFunctionSpace);

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = srcFieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  atlas::FieldSet tmpSrcFieldSet;
  atlas::FieldSet hybridFieldSet;
  for (auto & fieldname : sortedFieldNames) {
    if (variables.has(fieldname)) {
      atlas::Field hybridField =
        hybridFunctionSpace.createField<double>(
                              atlas::option::name(fieldname) |
                              atlas::option::levels(srcFieldSet[fieldname].shape(1)));
      atlas::Field tmpSrcField =
        srcFieldSet[fieldname].functionspace().createField<double>(
                                 atlas::option::name(fieldname) |
                                 atlas::option::levels(srcFieldSet[fieldname].shape(1)));
      atlas::array::make_view<double, 2>(hybridField).assign(0.0);

      const auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[fieldname]);
      auto srcTmpView = atlas::array::make_view<double, 2>(tmpSrcField);
      for (atlas::idx_t t = 0; t < srcFieldSet[fieldname].shape(0); ++t) {
        for (atlas::idx_t k = 0; k < srcFieldSet[fieldname].shape(1); ++k) {
          srcTmpView(t, k) = srcView(t, k);
        }
      }

      hybridFieldSet.add(hybridField);
      tmpSrcFieldSet.add(tmpSrcField);
    }
  }

  if (includingVectorInterpolation) {
    tmpSrcFieldSet["eastward_wind"].metadata().set("vector_field_name", "wind");
    tmpSrcFieldSet["northward_wind"].metadata().set("vector_field_name", "wind");
    hybridFieldSet["eastward_wind"].metadata().set("vector_field_name", "wind");
    hybridFieldSet["northward_wind"].metadata().set("vector_field_name", "wind");
  }


  interp.execute(tmpSrcFieldSet, hybridFieldSet);

  // Copy into StructuredColumns
  for (auto & fieldname : sortedFieldNames) {
    if (variables.has(fieldname)) {
      atlas::Field gaussField =
        gaussFunctionSpace.createField<double>(
                             atlas::option::name(fieldname) |
                             atlas::option::levels(srcFieldSet[fieldname].shape(1)));
      atlas::array::make_view<double, 2>(gaussField).assign(
        atlas::array::make_view<const double, 2>(hybridFieldSet[fieldname]));
      gaussField.set_dirty();  // atlas interpolation above produces dirty halos
      newFieldSet.add(gaussField);
    }
  }
}

// -----------------------------------------------------------------------------

Rescaling initRescaling(const GaussToCSParameters & params,
                        const eckit::mpi::Comm & comm,
                        const oops::Variables & activeVars,
                        const atlas::FunctionSpace & innerFspace,
                        const atlas::FunctionSpace & outerFspace,
                        const saber::interpolation::AtlasInterpWrapper & interp) {
  oops::Log::trace() << "::GaussToCS::initRescaling" << std::endl;
  if (params.interpolationRescaling.value().is_initialized()) {
    const auto & conf = params.interpolationRescaling.value().value();
    if (conf.has("horizontal covariance profile file path")) {
      // Compute rescaling fields from input covariance profile
      return Rescaling(comm,
                       conf,
                       activeVars,
                       innerFspace,
                       outerFspace,
                       interp);
    } else if (conf.has("input atlas file")) {
      // Read rescaling field from atlas file
      const auto readConf = conf.getSubConfiguration("input atlas file");
      return Rescaling(comm,
                       readConf,
                       activeVars,
                       outerFspace);
    } else {
      throw eckit::UserError("Missing parameters to initialize a relevant Rescaling object",
                             Here());
    }
  } else {
    return Rescaling();
  }
}
}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussToCS> makerGaussToCS_("gauss to cubed-sphere-dual");

// -----------------------------------------------------------------------------
// Note that this is slower than this needs to be
// as we are need to create 2 grid objects (very slow)
// In the future it might make sense to include an atlas grid (if available) from
// the model in outerGeometryData.
GaussToCS::GaussToCS(const oops::GeometryData & outerGeometryData,
                     const oops::Variables & outerVars,
                     const eckit::Configuration & covarConf,
                     const Parameters_ & params,
                     const oops::FieldSet3D & xb,
                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerVars_(outerVars),
    activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
    CSFunctionSpace_(outerGeometryData.functionSpace()),
    gaussGrid_(params.gaussGridUid.value()),
    interpType_(params.interpType.value()),
    gaussFunctionSpace_(createGaussFunctionSpace(gaussGrid_,
                                                 getHalo(params.interpType.value()))),
    gaussPartitioner_(new TransPartitioner()),
    csgrid_(CSFunctionSpace_.mesh().grid()),
    includingVectorInterpolation_(
         activeVars_.has("eastward_wind") && activeVars_.has("northward_wind")
         && params.skipVectorInterpolation == false ? true:false),
    interp_(gaussPartitioner_, gaussFunctionSpace_, csgrid_, CSFunctionSpace_,
            params.interpType.value(), includingVectorInterpolation_),
    inverseInterpolation_(createInverseInterpolation(
                              params.initializeInverseInterpolation.value(),
                              outerGeometryData.comm().size() == 1,
                              CSFunctionSpace_, gaussGrid_,
                              gaussPartitioner_,
                              interpType_,
                              includingVectorInterpolation_)),
    rescaling_(initRescaling(params,
                             outerGeometryData.comm(),
                             activeVars_,
                             gaussFunctionSpace_,
                             CSFunctionSpace_,
                             interp_)),

    innerGeometryData_(gaussFunctionSpace_, outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(),
                       outerGeometryData.comm())
{
  oops::Log::trace() << classname() << "::GaussToCS starting" << std::endl;
  oops::Log::trace() << classname() << "::GaussToCS done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = fieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // copy "passive variables"
  for (auto & fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      gaussFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // On input: fieldset on gaussian mesh

  // Create fieldset on cubed-sphere mesh.
  atlas::FieldSet csFieldSet;
  for (auto & fieldname : sortedFieldNames) {
      if (activeVars_.has(fieldname)) {
          atlas::Field csField =
              CSFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                    atlas::option::levels(gaussFieldSet[fieldname].shape(1)) |
                    atlas::option::halo(1));
          atlas::array::make_view<double, 2>(csField).assign(0.0);
          csFieldSet.add(csField);
      }
  }

  // Interpolate to cubed sphere
  interp_.execute(gaussFieldSet, csFieldSet);
  csFieldSet.set_dirty();  // atlas interpolation produces dirty halos


  for (auto & fieldname : sortedFieldNames) {
      if (activeVars_.has(fieldname)) {
        newFields.add(csFieldSet[fieldname]);
      }
  }

  fieldSet.fieldSet() = newFields;

  // Apply optional rescaling
  rescaling_.execute(fieldSet);

  oops::Log::trace() << classname() << "::multiply done"
                     << fieldSet.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // extract copy of vector of strings from fieldSet.field_names()
  // and apply sorting algorithm to copy
  auto sortedFieldNames = fieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // Apply optional rescaling
  rescaling_.execute(fieldSet);

  // Create empty fieldSets
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet csFieldSet = atlas::FieldSet();

  // copy "passive variables"
  for (auto & fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      csFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // Create gauss fieldset
  atlas::FieldSet gaussFieldSet;
  for (auto & fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      atlas::Field gaussField =
        gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
              atlas::option::levels(csFieldSet[fieldname].shape(1)) |
              atlas::option::halo(getHalo(interpType_)));
      atlas::array::make_view<double, 2>(gaussField).assign(0.0);
      gaussFieldSet.add(gaussField);
    }
  }

  // Adjoint of interpolation from gauss to dual cubed sphere
  interp_.executeAdjoint(gaussFieldSet, csFieldSet);

  for (auto & fieldname : sortedFieldNames) {
      if (activeVars_.has(fieldname)) {
          newFields.add(gaussFieldSet[fieldname]);
      }
  }

  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  atlas::FieldSet newFieldSet = atlas::FieldSet();
  atlas::FieldSet srcFieldSet = atlas::FieldSet();

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = fieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // copy "passive variables"
  for (auto & fieldname : sortedFieldNames) {
     if (activeVars_.has(fieldname)) {
       srcFieldSet.add(fieldSet[fieldname]);
     } else {
       newFieldSet.add(fieldSet[fieldname]);
     }
  }

  if (innerGeometryData_.comm().size() >= 2) {
    inverseInterpolateMultiplePEs(activeVars_, inverseInterpolation_,
                                  gaussFunctionSpace_,
                                  includingVectorInterpolation_,
                                  srcFieldSet, newFieldSet);
  } else {
    // A faster and more direct route is possible on a single PE
    inverseInterpolateSinglePE(activeVars_,
                               CSFunctionSpace_,
                               gaussFunctionSpace_,
                               gaussGrid_,
                               includingVectorInterpolation_,
                               interpType_,
                               srcFieldSet, newFieldSet);
  }

  fieldSet.fieldSet() = newFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::directCalibration(const oops::FieldSets & fset) {
  oops::Log::trace() << classname() << "::directCalibration start" << std::endl;

  oops::Log::trace() << classname() << "::directCalibration end" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCS::generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                                  const oops::Variables & innerVars) const {
  oops::FieldSet3D fset(this->validTime(), innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCS::generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                                  const oops::Variables & outerVars) const {
  oops::FieldSet3D fset(this->validTime(), outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void GaussToCS::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
