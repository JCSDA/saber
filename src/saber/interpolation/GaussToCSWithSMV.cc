/*
 * (C) Crown Copyright 2026- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/GaussToCSWithSMV.h"

#include <optional>

#include <algorithm>
#include <vector>

#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "eckit/exception/Exceptions.h"
#include "saber/interpolation/SMVInterpWrapper.h"


namespace saber {
namespace interpolation {

namespace {

int getHalo(const std::string& interpType) {
  oops::Log::trace() << "::GaussToCSWithSMV::getHalo" << std::endl;

  // Assumes that only interpolation types with "cubic" in their name require a
  // halo of 2.
  const std::string pattern = "cubic";
  if (interpType.find(pattern) != interpType.npos) {
    return 2;
  }
  return 1;
}

atlas::functionspace::StructuredColumns createGaussFunctionSpace(
    const atlas::StructuredGrid& gaussGrid,
    const atlas::grid::Partitioner& partitioner, const int halo) {
  oops::Log::trace() << "::GaussToCSWithSMV::createGaussFunctionSpace" << std::endl;

  const auto gaussFSpace = atlas::functionspace::StructuredColumns(
      gaussGrid, partitioner, atlas::option::halo(halo));

  return gaussFSpace;
}

// -----------------------------------------------------------------------------

Rescaling initRescaling(
    const GaussToCSWithSMVParameters& params, const eckit::mpi::Comm& comm,
    const oops::Variables& activeVars, const atlas::FunctionSpace& innerFspace,
    const atlas::FunctionSpace& outerFspace,
    const saber::interpolation::AtlasInterpWrapper& interp) {
  oops::Log::trace() << "::GaussToCS::initRescaling" << std::endl;
  if (params.interpolationRescaling.value().is_initialized()) {
    const auto& conf = params.interpolationRescaling.value().value();
    if (conf.has("horizontal covariance profile file path")) {
      // Compute rescaling fields from input covariance profile
      return Rescaling(comm, conf, activeVars, innerFspace, outerFspace,
                       interp);
    } else if (conf.has("input atlas file")) {
      // Read rescaling field from atlas file
      const auto readConf = conf.getSubConfiguration("input atlas file");
      return Rescaling(comm, readConf, activeVars, outerFspace);
    } else {
      throw eckit::UserError(
          "Missing parameters to initialize a relevant Rescaling object",
          Here());
    }
  } else {
    return Rescaling();
  }
}

// Util for gathering a subset of field configurations, when creating fields for
// interpolation.
std::vector<atlas::util::Config> gatherInterpFieldConfigs(
    const atlas::FieldSet& fset, const oops::Variables& variables) {
  std::vector<atlas::util::Config> out;
  for (auto& fieldname : fset.field_names()) {
    if (variables.has(fieldname)) {
      out.emplace_back(atlas::option::name(fieldname) |
                       atlas::option::levels(fset[fieldname].shape(1)));
    }
  }
  return out;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussToCSWithSMV> makerGaussToCSWithSMV_(
    "gauss to cubed-sphere-dual with smv-interp");

// -----------------------------------------------------------------------------
// Note that this is slower than this needs to be
// as we are need to create 2 grid objects (very slow)
// In the future it might make sense to include an atlas grid (if available)
// from the model in outerGeometryData.
GaussToCSWithSMV::GaussToCSWithSMV(const oops::GeometryData& outerGeometryData,
                                   const oops::Variables& outerVars,
                                   const eckit::Configuration& covarConf,
                                   const Parameters_& params,
                                   const oops::FieldSet3D& xb,
                                   const oops::FieldSet3D& fg)
    : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
      innerVars_(outerVars),
      activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
      CSFunctionSpace_(outerGeometryData.functionSpace()),
      gaussGrid_(params.gaussGridUid.value()),
      interpType_(params.interpType.value()),
      gaussPartitioner_([&]() {
        eckit::mpi::setCommDefault(outerGeometryData.comm().name());
        return new atlas::grid::detail::partitioner::TransPartitioner();
      }()),
      gaussFunctionSpace_(createGaussFunctionSpace(
          gaussGrid_, gaussPartitioner_, getHalo(params.interpType.value()))),
      csgrid_(CSFunctionSpace_.mesh().grid()),
      includingVectorInterpolation_(params.interpolateWindAsVector.value()),
      interp_([&]() {
        eckit::mpi::setCommDefault(outerGeometryData.comm().name());
        return saber::interpolation::AtlasInterpWrapper(
            gaussPartitioner_, gaussFunctionSpace_, csgrid_, CSFunctionSpace_,
            params.interpType.value(), includingVectorInterpolation_);
      }()),
      inverseInterpolation_([&]() -> std::optional<SMVInterpWrapper> {
        if (params.initializeInverseInterpolation.value()) {
          return SMVInterpWrapper(CSFunctionSpace_,
                                  gaussFunctionSpace_,
                                  gaussGrid_,
                                  includingVectorInterpolation_,
                                  params.normalizeInverseInterpolation.value());
        } else {
          return std::nullopt;
        }
      }()),
      rescaling_(initRescaling(params, outerGeometryData.comm(), activeVars_,
                               gaussFunctionSpace_, CSFunctionSpace_, interp_)),

      innerGeometryData_(gaussFunctionSpace_, outerGeometryData.fieldSet(),
                         outerGeometryData.levelsAreTopDown(),
                         outerGeometryData.comm()) {
  oops::Log::trace() << classname() << "::GaussToCS starting" << std::endl;
  oops::Log::trace() << classname() << "::GaussToCS done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCSWithSMV::multiply(oops::FieldSet3D& fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = fieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // copy "passive variables"
  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      gaussFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // On input: fieldset on gaussian mesh

  // Create fieldset on cubed-sphere mesh.
  atlas::FieldSet csFieldSet;
  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      atlas::Field csField = CSFunctionSpace_.createField<double>(
          atlas::option::name(fieldname) |
          atlas::option::levels(gaussFieldSet[fieldname].shape(1)) |
          atlas::option::halo(1));
      atlas::array::make_view<double, 2>(csField).assign(0.0);
      csFieldSet.add(csField);
    }
  }

  // Interpolate to cubed sphere
  interp_.execute(gaussFieldSet, csFieldSet);
  csFieldSet.set_dirty();  // atlas interpolation produces dirty halos

  for (auto& fieldname : sortedFieldNames) {
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

void GaussToCSWithSMV::multiplyAD(oops::FieldSet3D& fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

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
  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      csFieldSet.add(fieldSet[fieldname]);
    } else {
      newFields.add(fieldSet[fieldname]);
    }
  }

  // Create gauss fieldset
  atlas::FieldSet gaussFieldSet;
  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      atlas::Field gaussField = gaussFunctionSpace_.createField<double>(
          atlas::option::name(fieldname) |
          atlas::option::levels(csFieldSet[fieldname].shape(1)) |
          atlas::option::halo(getHalo(interpType_)));
      atlas::array::make_view<double, 2>(gaussField).assign(0.0);
      gaussFieldSet.add(gaussField);
    }
  }

  // Adjoint of interpolation from gauss to dual cubed sphere
  interp_.executeAdjoint(gaussFieldSet, csFieldSet);

  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      newFields.add(gaussFieldSet[fieldname]);
    }
  }
  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCSWithSMV::inverseMultiply(oops::FieldSet3D& fieldSet) const {
  if (!inverseInterpolation_) {
    std::string errMsg =
        "inverseInterpolation_ optional not initialised to CS2GaussWithSMV. "
        "Aborting.";
    throw eckit::Exception(errMsg, Here());
  }

  oops::Log::trace() << classname() << "::inverseMultiply starting"
                     << std::endl;

  atlas::FieldSet newFieldSet = atlas::FieldSet();
  atlas::FieldSet srcFieldSet = atlas::FieldSet();

  // extract copy of field names and apply sorting algorithm
  auto sortedFieldNames = fieldSet.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // copy "passive variables"
  for (auto& fieldname : sortedFieldNames) {
    if (activeVars_.has(fieldname)) {
      srcFieldSet.add(fieldSet[fieldname]);
    } else {
      newFieldSet.add(fieldSet[fieldname]);
    }
  }

  const std::vector<atlas::util::Config> fieldConfigs =
      gatherInterpFieldConfigs(srcFieldSet, activeVars_);
  inverseInterpolation_->doInterpolation(fieldConfigs, srcFieldSet, newFieldSet);

  fieldSet.fieldSet() = newFieldSet;

  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCSWithSMV::directCalibration(const oops::FieldSets& fset) {
  oops::Log::trace() << classname() << "::directCalibration start" << std::endl;

  oops::Log::trace() << classname() << "::directCalibration end" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCSWithSMV::generateInnerFieldSet(
    const oops::GeometryData& innerGeometryData,
    const oops::Variables& innerVars) const {
  oops::FieldSet3D fset(validTime_, innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(
      innerGeometryData.comm(), innerGeometryData.functionSpace(), innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D GaussToCSWithSMV::generateOuterFieldSet(
    const oops::GeometryData& outerGeometryData,
    const oops::Variables& outerVars) const {
  oops::FieldSet3D fset(validTime_, outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(
      outerGeometryData.comm(), outerGeometryData.functionSpace(), outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void GaussToCSWithSMV::print(std::ostream& os) const { os << classname(); }

}  // namespace interpolation
}  // namespace saber
