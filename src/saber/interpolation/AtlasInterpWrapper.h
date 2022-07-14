/*
 * (C) Crown Copyright 2021-2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_INTERPOLATION_ATLASINTERPWRAPPER_H_
#define SABER_INTERPOLATION_ATLASINTERPWRAPPER_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/redistribution/Redistribution.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace atlas {
class Field;
class FieldSet;
}

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

class AtlasInterpWrapper {

 public:
  static const std::string classname() {return "saber::interpolation::AtlasInterpWrapper";}

  AtlasInterpWrapper(const atlas::grid::Partitioner &,
                     const atlas::FunctionSpace &,
                     const atlas::Grid &,
                     const atlas::FunctionSpace &);
  AtlasInterpWrapper() {}

  void execute(const atlas::Field & srcField,
               atlas::Field & targetField) const;

  void executeAdjoint(atlas::Field & srcField,
                      const atlas::Field & targetField) const;

  void execute(const atlas::FieldSet & srcFieldSet,
               atlas::FieldSet & targetFieldSet) const;

  void executeAdjoint(atlas::FieldSet & srcFieldSet,
                      const atlas::FieldSet & targetFieldSet) const;

 private:
  void print(std::ostream &) const {}
  atlas::FunctionSpace targetFunctionSpace_; 
  atlas::Interpolation interp_;
  atlas::Redistribution redistr_;
};

}  // namespace interpolation
}  // namespace saber

// -----------------------------------------------------------------------------

namespace saber {
namespace interpolation {
namespace detail {

// -----------------------------------------------------------------------------

atlas::FunctionSpace createTargetFunctionSpace(
  const atlas::grid::Partitioner & sourcePartitioner,
  const atlas::FunctionSpace & sourceFunctionSpace,
  const atlas::Grid & outputGrid,
  const std::string & outputFunctionSpaceType) {
  oops::Log::trace() << "createTargetFunctionSpace starting" << std::endl;

  // Get or compute source mesh
  atlas::Mesh sourceMesh;
  if (sourceFunctionSpace.type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(sourceFunctionSpace);
    sourceMesh = atlas::MeshGenerator("structured").generate(fs.grid(), sourcePartitioner);
  } else if (sourceFunctionSpace.type() == "NodeColumns") {
    atlas::functionspace::NodeColumns fs(sourceFunctionSpace);
    sourceMesh = fs.mesh();
  } else {
    ABORT(sourceFunctionSpace.type() + " source function space not supported yet");
  }

  // Get target partitioner from source mesh
  atlas::grid::MatchingMeshPartitioner targetPartitioner(sourceMesh);

  // Create target function space
  if (outputFunctionSpaceType == "StructuredColumns") {
    // StructuredColumns
    const atlas::grid::Distribution distribution(outputGrid, targetPartitioner);
    return atlas::functionspace::StructuredColumns(outputGrid, distribution);
  } else if (outputFunctionSpaceType == "NodeColumns") {
    // NodeColumns
    if (outputGrid.name().substr(0, 2).compare("CS") == 0) {
      // CubedSphere
      atlas::Mesh targetMesh = atlas::MeshGenerator("cubedsphere").generate(outputGrid,
        targetPartitioner);
      return atlas::functionspace::CubedSphereNodeColumns(targetMesh);
    } else {
      // Other NodeColumns (TODO: not working!)
      atlas::Mesh targetMesh = atlas::MeshGenerator("delaunay").generate(outputGrid,
        targetPartitioner);
      return atlas::functionspace::NodeColumns(targetMesh);
    }
  } else {
    ABORT(outputFunctionSpaceType + " output function space type not supported yet");
    return atlas::FunctionSpace();
  }

  oops::Log::trace() << "createTargetFunctionSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

atlas::Interpolation createAtlasInterpolation(const atlas::FunctionSpace & inputFS,
                                              const atlas::FunctionSpace & matchingFS) {
  oops::Log::trace() << "createAtlasInterpolation starting" << std::endl;

  atlas::util::Config interpConfig;
  if (inputFS.type() == "StructuredColumns") {
    // StructuredColumns
    interpConfig.set("type", "structured-linear2D");
  } else if (inputFS.type() == "NodeColumns") {
    // NodeColumns
    interpConfig.set("type", "unstructured-bilinear-lonlat");
  } else {
    ABORT(inputFS.type() + " source function space type not supported yet");
  }
  interpConfig.set("adjoint", "true");

  atlas::Interpolation interp(interpConfig, inputFS, matchingFS);
  oops::Log::trace() << "createAtlasInterpolation done" << std::endl;
  return interp;
}

// -----------------------------------------------------------------------------

atlas::Redistribution createAtlasRedistribution(
    const atlas::FunctionSpace & matchingFS,
    const atlas::FunctionSpace & outputFS) {
  oops::Log::trace() << "createAtlasRedistribution starting" << std::endl;
  atlas::Redistribution redistr(matchingFS, outputFS);
  oops::Log::trace() << "createAtlasRedistribution done" << std::endl;
  return redistr;
}

// -----------------------------------------------------------------------------

}  // namespace detail
}  // namespace interpolation
}  // namespace saber

// ------------------------------------------------------------------------------------------------

namespace saber {
namespace interpolation {

AtlasInterpWrapper::AtlasInterpWrapper(const atlas::grid::Partitioner & sourcePartitioner,
                                       const atlas::FunctionSpace & sourceFunctionSpace,
                                       const atlas::Grid & outputGrid,
                                       const atlas::FunctionSpace & outputFunctionSpace) :
  targetFunctionSpace_(detail::createTargetFunctionSpace(sourcePartitioner, sourceFunctionSpace,
    outputGrid, outputFunctionSpace.type())),
  interp_(detail::createAtlasInterpolation(sourceFunctionSpace, targetFunctionSpace_)),
  redistr_(detail::createAtlasRedistribution(targetFunctionSpace_, outputFunctionSpace))
{oops::Log::trace() << "AtlasInterpWrapper::AtlasInterpWrapper done" << std::endl;};

//-------------------------------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::Field & srcField,
                                 atlas::Field & outputField) const {

  // Copy source field to exchange halo

  // Create empty source field
  atlas::Field srcTmp = srcField.functionspace().createField<double>
      (atlas::option::name(srcField.name()) |
       atlas::option::levels(srcField.levels()));

  // Copy source field
  auto src_v = atlas::array::make_view<double, 2>(srcField);
  auto srcTmp_v = atlas::array::make_view<double, 2>(srcTmp);
  for (atlas::idx_t t = 0; t < srcField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < srcField.shape(1); ++k) {
      srcTmp_v(t, k) = src_v(t, k);
    }
  }

  // Exchange halo
  srcTmp.haloExchange();

  // Create empty target field
  auto tmpField =
    targetFunctionSpace_.createField<double>(atlas::option::name(outputField.name()) |
                                             atlas::option::levels(srcField.levels()));

  // Initialize target field to zero
  auto tmp_v = atlas::array::make_view<double, 2>(tmpField);
  tmp_v.assign(0.0);

  // Interpolate from source field to target field
  interp_.execute(srcTmp, tmpField);

  // Redistribute from target field to output field
  redistr_.execute(tmpField, outputField);

  // Exchange halo
  outputField.haloExchange();
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::Field & srcField,
                                        const atlas::Field & outputField) const {
  oops::Log::trace() << "AtlasInterpWrapper::executeAdjoint start "
                        "srcFieldName outputFieldName" <<
                        srcField.name() << " " << outputField.name() << std::endl;

  // create a field that is a copy of target
  // take halo adjoint of temp field
  atlas::Field tmpField =
    outputField.functionspace().createField<double>(atlas::option::name(outputField.name()) |
                                                    atlas::option::levels(outputField.levels()));

  auto tmp_v = atlas::array::make_view<double, 2>(tmpField);
  auto tgt_v = atlas::array::make_view<double, 2>(outputField);
  for (atlas::idx_t t = 0; t < outputField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < outputField.shape(1); ++k) {
      tmp_v(t, k) = tgt_v(t, k);
    }
  }
  tmpField.adjointHaloExchange();

  auto tmp2Field = targetFunctionSpace_.createField<double>(
                          atlas::option::name(outputField.name()) |
                          atlas::option::levels(outputField.levels()));
  auto tmp2_v = atlas::array::make_view<double, 2>(tmp2Field);
  tmp2_v.assign(0.0);

  redistr_.execute(tmpField, tmp2Field);

  auto src_v = atlas::array::make_view<double, 2>(srcField);
  src_v.assign(0.0);

  interp_.execute_adjoint(srcField, tmp2Field);

  oops::Log::trace() << "AtlasInterpWrapper::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::FieldSet & srcFieldSet,
                                 atlas::FieldSet & targetFieldSet) const {
  for (auto & srcField : srcFieldSet) {
    execute(srcField, targetFieldSet[srcField.name()]);
  }
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                               const atlas::FieldSet & targetFieldSet) const {
  for (auto & srcField : srcFieldSet) {
    executeAdjoint(srcField, targetFieldSet[srcField.name()]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace saber

// ------------------------------------------------------------------------------------------------

#endif  // SABER_INTERPOLATION_ATLASINTERPWRAPPER_H_
