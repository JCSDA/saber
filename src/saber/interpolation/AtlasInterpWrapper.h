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
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Point.h"
#include "atlas/util/PolygonLocator.h"
#include "atlas/util/PolygonXY.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"


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
  ~AtlasInterpWrapper() {}

  void execute(const atlas::Field &, atlas::Field &) const;

  void executeAdjoint(atlas::Field &, const atlas::Field &) const;

  void execute(const atlas::FieldSet &, atlas::FieldSet &) const;

  void executeAdjoint(atlas::FieldSet &, const atlas::FieldSet &) const;

 private:
  void print(std::ostream &) const {}
  atlas::FunctionSpace localOutputFunctionSpace_;
  std::vector<size_t> localTask_;
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
  const atlas::grid::Partitioner & srcPartitioner,
  const atlas::FunctionSpace & srcFunctionSpace,
  const atlas::Grid & outputGrid,
  const std::string & outputFunctionSpaceType) {
  oops::Log::trace() << "createTargetFunctionSpace starting" << std::endl;

  // Check source partitioner type
  if (srcPartitioner.type() != "trans") {
    ABORT("createTargetFunctionSpace: trans partitioner required");
  }

  // Get or compute source mesh
  atlas::Mesh srcMesh;
  if (srcFunctionSpace.type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(srcFunctionSpace);
    srcMesh = atlas::MeshGenerator("structured").generate(fs.grid(), srcPartitioner);
  } else if (srcFunctionSpace.type() == "NodeColumns") {
    atlas::functionspace::NodeColumns fs(srcFunctionSpace);
    srcMesh = fs.mesh();
  } else {
    ABORT(srcFunctionSpace.type() + " source function space not supported yet");
  }

  // Get target partitioner from source mesh
  atlas::grid::MatchingMeshPartitioner targetPartitioner(srcMesh);

  // Create target function space
  atlas::FunctionSpace targetFunctionSpace;
  if (outputFunctionSpaceType == "StructuredColumns") {
    // StructuredColumns
    const atlas::grid::Distribution distribution(outputGrid, targetPartitioner);
    targetFunctionSpace = atlas::functionspace::StructuredColumns(outputGrid, distribution);
  } else if (outputFunctionSpaceType == "NodeColumns") {
    // NodeColumns
    if (outputGrid.name().substr(0, 2).compare("CS") == 0) {
      // CubedSphere
      atlas::Mesh targetMesh = atlas::MeshGenerator("cubedsphere").generate(outputGrid,
        targetPartitioner);
      targetFunctionSpace = atlas::functionspace::CubedSphereNodeColumns(targetMesh);
    } else {
      // Other NodeColumns (TODO: not working!)
      atlas::Mesh targetMesh = atlas::MeshGenerator("delaunay").generate(outputGrid, targetPartitioner);
      targetFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);
    }
  } else {
    ABORT(outputFunctionSpaceType + " output function space type not supported yet");
    targetFunctionSpace = atlas::FunctionSpace();
  }

  oops::Log::trace() << "createTargetFunctionSpace done" << std::endl;
  return targetFunctionSpace;
}

// -----------------------------------------------------------------------------

atlas::Interpolation createAtlasInterpolation(const atlas::FunctionSpace & inputFunctionSpace,
                                              const atlas::FunctionSpace & matchingFunctionSpace) {
  oops::Log::trace() << "createAtlasInterpolation starting" << std::endl;

  atlas::util::Config interpConfig;
  if (inputFunctionSpace.type() == "StructuredColumns") {
    // StructuredColumns
    interpConfig.set("type", "structured-linear2D");
  } else if (inputFunctionSpace.type() == "NodeColumns") {
    // NodeColumns
    interpConfig.set("type", "unstructured-bilinear-lonlat");
  } else {
    ABORT(inputFunctionSpace.type() + " source function space type not supported yet");
  }
  interpConfig.set("adjoint", "true");

  atlas::Interpolation interp(interpConfig, inputFunctionSpace, matchingFunctionSpace);
  oops::Log::trace() << "createAtlasInterpolation done" << std::endl;
  return interp;
}

// -----------------------------------------------------------------------------

atlas::Redistribution createAtlasRedistribution(
    const atlas::FunctionSpace & matchingFunctionSpace,
    const atlas::FunctionSpace & outputFunctionSpace) {
  oops::Log::trace() << "createAtlasRedistribution starting" << std::endl;
  atlas::Redistribution redistr(matchingFunctionSpace, outputFunctionSpace);
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

AtlasInterpWrapper::AtlasInterpWrapper(const atlas::grid::Partitioner & srcPartitioner,
                                       const atlas::FunctionSpace & srcFunctionSpace,
                                       const atlas::Grid & outputGrid,
                                       const atlas::FunctionSpace & outputFunctionSpace) :
  localOutputFunctionSpace_(), targetFunctionSpace_(), interp_(), redistr_() {

  oops::Log::trace() << "AtlasInterpWrapper::AtlasInterpWrapper starting" << std::endl;

  if (outputFunctionSpace.type() == "PointCloud") {
    // PointCloud output function space
    atlas::functionspace::PointCloud fs(outputFunctionSpace);

    // Split points between distribution polygons
    atlas::util::ListPolygonXY polygons(srcFunctionSpace.polygons());
    atlas::util::PolygonLocator find_partition(polygons);
    std::vector<atlas::PointXY> points;
    auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
    for (atlas::idx_t i = 0; i < fs.size(); ++i) {
      atlas::PointLonLat pointLonLat(lonlatView(i, 0), lonlatView(i, 1));
      pointLonLat.normalise();
      atlas::PointXY point(pointLonLat);
      localTask_.push_back(find_partition(point));
      if (eckit::mpi::comm().rank() == localTask_[i]) {
        points.push_back(point);
      }
    }
    localOutputFunctionSpace_ = atlas::functionspace::PointCloud(points);   

    // Create interpolation
    if (localOutputFunctionSpace_.size() > 0) {
      interp_ = detail::createAtlasInterpolation(srcFunctionSpace, localOutputFunctionSpace_);
    }
  } else {
    // Other output function space

    // Create target function space
    targetFunctionSpace_ = detail::createTargetFunctionSpace(srcPartitioner, srcFunctionSpace,
      outputGrid, outputFunctionSpace.type());

    // Create interpolation
    interp_ = detail::createAtlasInterpolation(srcFunctionSpace, targetFunctionSpace_);

    // Create redistribution
    redistr_ = detail::createAtlasRedistribution(targetFunctionSpace_, outputFunctionSpace);
  }

  oops::Log::trace() << "AtlasInterpWrapper::AtlasInterpWrapper done" << std::endl;
}

//-------------------------------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::Field & srcField,
                                 atlas::Field & outputField) const {

  // Copy source field to exchange halo

  // Create empty source field
  atlas::Field srcTmp = srcField.functionspace().createField<double>(
    atlas::option::name(srcField.name()) | atlas::option::levels(srcField.levels()));

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

  if (outputField.functionspace().type() == "PointCloud") {
    // PointCloud output function space
    atlas::functionspace::PointCloud fs(outputField.functionspace());
    atlas::functionspace::PointCloud localFs(localOutputFunctionSpace_);

    // Initialize global vector
    std::vector<double> globalData(fs.size()*outputField.levels(), 0.0);

    if (localFs.size() > 0) {
      // Create local output field
      atlas::Field localOutputField = localOutputFunctionSpace_.createField<double>(
      atlas::option::name(outputField.name()) | atlas::option::levels(outputField.levels()));

      // Interpolate from source field to local output field
      interp_.execute(srcTmp, localOutputField);

      // Copy local output field into global vector
      auto localOutputView = atlas::array::make_view<double, 2>(localOutputField);
      atlas::idx_t globalInc = 0;
      for (atlas::idx_t k = 0; k < outputField.levels(); ++k) {
        atlas::idx_t localInc = 0;
        for (atlas::idx_t i = 0; i < fs.size(); ++i) {
          if (eckit::mpi::comm().rank() == localTask_[i]) {
            // Copy local value
            globalData[globalInc] = localOutputView(localInc, k);
            localInc += 1;
          }
          globalInc += 1;
        }
      }
    }

    // Sum global vector over MPI tasks
    eckit::mpi::comm().allReduceInPlace(globalData.begin(), globalData.end(), eckit::mpi::sum());

    // Copy global vector into global field
    auto outputView = atlas::array::make_view<double, 2>(outputField);
    atlas::idx_t globalInc = 0;
    for (atlas::idx_t k = 0; k < outputField.levels(); ++k) {
      for (atlas::idx_t i = 0; i < fs.size(); ++i) {
        outputView(i, k) = globalData[globalInc];
        globalInc += 1;
      }
    }
  } else {
    // Other output function space

    // Create empty target field
    auto targetField =
      targetFunctionSpace_.createField<double>(atlas::option::name(outputField.name()) |
                                            atlas::option::levels(srcField.levels()));

    // Initialize target field to zero
    auto tmp_v = atlas::array::make_view<double, 2>(targetField);
    tmp_v.assign(0.0);

    // Interpolate from source field to target field
    interp_.execute(srcTmp, targetField);

    // Redistribute from target field to output field
    redistr_.execute(targetField, outputField);

    // Exchange halo
    outputField.haloExchange();
  }
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::Field & srcField,
                                        const atlas::Field & outputField) const {
  oops::Log::trace() << "AtlasInterpWrapper::executeAdjoint start "
                        "srcFieldName outputFieldName" <<
                        srcField.name() << " " << outputField.name() << std::endl;

  // create a field that is a copy of target
  // take halo adjoint of temp field
  atlas::Field targetField = outputField.functionspace().createField<double>(
    atlas::option::name(outputField.name()) | atlas::option::levels(outputField.levels()));

  auto tmp_v = atlas::array::make_view<double, 2>(targetField);
  auto target_v = atlas::array::make_view<double, 2>(outputField);
  for (atlas::idx_t t = 0; t < outputField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < outputField.shape(1); ++k) {
      tmp_v(t, k) = target_v(t, k);
    }
  }
  targetField.adjointHaloExchange();

  auto tmp2Field = targetFunctionSpace_.createField<double>(
                          atlas::option::name(outputField.name()) |
                          atlas::option::levels(outputField.levels()));
  auto tmp2_v = atlas::array::make_view<double, 2>(tmp2Field);
  tmp2_v.assign(0.0);

  redistr_.execute(targetField, tmp2Field);

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
