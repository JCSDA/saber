/*
 * (C) Crown Copyright 2021-2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/AtlasInterpWrapper.h"

#include "eckit/exception/Exceptions.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace interpolation {
namespace  {

// -----------------------------------------------------------------------------

atlas::FunctionSpace createTargetFunctionSpace(
  const atlas::grid::Partitioner & srcPartitioner,
  const atlas::FunctionSpace & srcFunctionSpace,
  const atlas::Grid & dstGrid,
  const std::string & dstFunctionSpaceType) {
  oops::Log::trace() << "createTargetFunctionSpace starting" << std::endl;

  // Get or compute source mesh
  atlas::Mesh srcMesh;
  if (srcFunctionSpace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns fs(srcFunctionSpace);
    srcMesh = atlas::MeshGenerator("structured").generate(fs.grid(), srcPartitioner);
  } else if (srcFunctionSpace.type() == "NodeColumns") {
    const atlas::functionspace::NodeColumns fs(srcFunctionSpace);
    srcMesh = fs.mesh();
  } else {
    throw eckit::FunctionalityNotSupported(srcFunctionSpace.type()
      + " source function space not supported yet", Here());
  }

  // Get target partitioner from source mesh
  const atlas::grid::MatchingMeshPartitioner targetPartitioner(srcMesh);

  // Create target function space
  atlas::FunctionSpace targetFunctionSpace;
  if (dstFunctionSpaceType == "StructuredColumns") {
    // StructuredColumns
    const atlas::grid::Distribution distribution(dstGrid, targetPartitioner);
    targetFunctionSpace = atlas::functionspace::StructuredColumns(dstGrid, distribution);
  } else if (dstFunctionSpaceType == "NodeColumns") {
    // NodeColumns
    if (dstGrid.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      const atlas::Mesh targetMesh = atlas::MeshGenerator("cubedsphere_dual").generate(dstGrid,
        targetPartitioner);
      targetFunctionSpace = atlas::functionspace::CubedSphereNodeColumns(targetMesh);
    } else {
      // Other NodeColumns (TODO: not working!)
      const atlas::Mesh targetMesh = atlas::MeshGenerator("delaunay").generate(dstGrid,
        targetPartitioner);
      targetFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);
    }
  } else {
    throw eckit::FunctionalityNotSupported(dstFunctionSpaceType
      + " destination function space type not supported yet", Here());
    targetFunctionSpace = atlas::FunctionSpace();
  }

  oops::Log::trace() << "createTargetFunctionSpace done" << std::endl;
  return targetFunctionSpace;
}

// -----------------------------------------------------------------------------

atlas::Interpolation createAtlasInterpolation(const atlas::FunctionSpace & srcFunctionSpace,
                                              const atlas::FunctionSpace & matchingFunctionSpace) {
  oops::Log::trace() << "createAtlasInterpolation starting" << std::endl;

  atlas::util::Config interpConfig;
  if (srcFunctionSpace.type() == "StructuredColumns") {
    // StructuredColumns
    interpConfig.set("type", "structured-linear2D");
  } else if (srcFunctionSpace.type() == "NodeColumns") {
    // NodeColumns
    interpConfig.set("type", "unstructured-bilinear-lonlat");
  } else {
    throw eckit::FunctionalityNotSupported(srcFunctionSpace.type()
      + " source function space type not supported yet", Here());
  }
  interpConfig.set("adjoint", "true");

  const atlas::Interpolation interp(interpConfig, srcFunctionSpace, matchingFunctionSpace);
  oops::Log::trace() << "createAtlasInterpolation done" << std::endl;
  return interp;
}

// -----------------------------------------------------------------------------

atlas::Redistribution createAtlasRedistribution(
    const atlas::FunctionSpace & matchingFunctionSpace,
    const atlas::FunctionSpace & dstFunctionSpace) {
  oops::Log::trace() << "createAtlasRedistribution starting" << std::endl;
  const atlas::Redistribution redistr(matchingFunctionSpace, dstFunctionSpace);
  oops::Log::trace() << "createAtlasRedistribution done" << std::endl;
  return redistr;
}

// -----------------------------------------------------------------------------

atlas::Redistribution createInverseAtlasRedistribution(
    const atlas::FunctionSpace & matchingFS,
    const atlas::FunctionSpace & outputFS) {

  return atlas::Redistribution(outputFS, matchingFS);
}

// -----------------------------------------------------------------------------

}  // namespace
}  // namespace interpolation
}  // namespace saber

// ------------------------------------------------------------------------------------------------

namespace saber {
namespace interpolation {

AtlasInterpWrapper::AtlasInterpWrapper(const atlas::grid::Partitioner & srcPartitioner,
                                       const atlas::FunctionSpace & srcFunctionSpace,
                                       const atlas::Grid & dstGrid,
                                       const atlas::FunctionSpace & dstFunctionSpace) :
  localDstFunctionSpace_(), targetFunctionSpace_(), interp_(), redistr_(), inverseRedistr_() {
  oops::Log::trace() << "AtlasInterpWrapper::AtlasInterpWrapper starting" << std::endl;

  if (dstFunctionSpace.type() == "PointCloud") {
    // PointCloud destination function space
    atlas::functionspace::PointCloud fs(dstFunctionSpace);

    // Splitting points between distribution polygons
    const atlas::util::ListPolygonXY polygons(srcFunctionSpace.polygons());
    const atlas::util::PolygonLocator find_partition(polygons);
    std::vector<atlas::PointXY> points;
    auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
    for (atlas::idx_t i = 0; i < fs.size(); ++i) {
      atlas::PointLonLat pointLonLat(lonlatView(i, 0), lonlatView(i, 1));
      pointLonLat.normalise();
      const atlas::PointXY point(pointLonLat);
      localTask_.push_back(find_partition(point));
      if (eckit::mpi::comm().rank() == localTask_[i]) {
        points.push_back(point);
      }
    }

    // Add extra points if needed
    auto ghostView = atlas::array::make_view<int, 1>(srcFunctionSpace.ghost());
    lonlatView = atlas::array::make_view<double, 2>(srcFunctionSpace.lonlat());
    for (atlas::idx_t jnode = 0; jnode < srcFunctionSpace.size(); ++jnode) {
      if (ghostView(jnode) == 0) {
        if (points.size() == 0) {
           atlas::PointLonLat pointLonLat(lonlatView(jnode, 0), lonlatView(jnode, 1));
           pointLonLat.normalise();
           const atlas::PointXY point(pointLonLat);
           points.push_back(point);
        }
      }
    }

    // Create local destination function space
    localDstFunctionSpace_ = atlas::functionspace::PointCloud(points);

    // Interpolation setup
    if (localDstFunctionSpace_.size() > 0) {
      interp_ = createAtlasInterpolation(srcFunctionSpace, localDstFunctionSpace_);
    }

  } else {
    // Other destination function space

    // Target function space setup
    targetFunctionSpace_ = createTargetFunctionSpace(srcPartitioner, srcFunctionSpace,
      dstGrid, dstFunctionSpace.type());

    // Interpolation
    interp_ = createAtlasInterpolation(srcFunctionSpace, targetFunctionSpace_);

    // Redistribution
    redistr_ = createAtlasRedistribution(targetFunctionSpace_, dstFunctionSpace);

    // Inverse Redistribution
    inverseRedistr_ = createInverseAtlasRedistribution(targetFunctionSpace_,
                                                               dstFunctionSpace);
  }

  oops::Log::trace() << "AtlasInterpWrapper::AtlasInterpWrapper done" << std::endl;
}

//-------------------------------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::Field & srcField,
                                 atlas::Field & dstField) const {
  // Empty source field setup
  atlas::Field srcTmpField = srcField.functionspace().createField(srcField);

  // Copy of source field
  const auto srcView = atlas::array::make_view<double, 2>(srcField);
  auto srcTmpView = atlas::array::make_view<double, 2>(srcTmpField);
  for (atlas::idx_t t = 0; t < srcField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < srcField.shape(1); ++k) {
      srcTmpView(t, k) = srcView(t, k);
    }
  }

  if (dstField.functionspace().type() == "PointCloud") {
    // PointCloud destination function space
    const atlas::functionspace::PointCloud fs(dstField.functionspace());
    const atlas::functionspace::PointCloud localFs(localDstFunctionSpace_);

    // Global vector initialization
    std::vector<double> globalData(fs.size()*dstField.levels(), 0.0);

    // Local destination field setup
    atlas::Field localDstField = localDstFunctionSpace_.createField<double>(
      atlas::option::name(dstField.name()) | atlas::option::levels(dstField.levels()));

    // Interpolation from source field to local destination field
    interp_.execute(srcTmpField, localDstField);

    // Copy of local destination field into global vector
    const auto localDstView = atlas::array::make_view<double, 2>(localDstField);
    atlas::idx_t globalInc = 0;
    for (atlas::idx_t k = 0; k < dstField.levels(); ++k) {
      atlas::idx_t localInc = 0;
      for (atlas::idx_t i = 0; i < fs.size(); ++i) {
        if (eckit::mpi::comm().rank() == localTask_[i]) {
          globalData[globalInc] = localDstView(localInc, k);
          localInc += 1;
        }
        globalInc += 1;
      }
    }


    // Sum of global vector over MPI tasks
    eckit::mpi::comm().allReduceInPlace(globalData.begin(), globalData.end(), eckit::mpi::sum());

    // Copy of global vector into global field
    auto dstView = atlas::array::make_view<double, 2>(dstField);
    globalInc = 0;
    for (atlas::idx_t k = 0; k < dstField.levels(); ++k) {
      for (atlas::idx_t i = 0; i < fs.size(); ++i) {
        dstView(i, k) = globalData[globalInc];
        globalInc += 1;
      }
    }
  } else {
    // Other destination function space

    // Empty target field setup
    auto targetField =
      targetFunctionSpace_.createField<double>(atlas::option::name(dstField.name()) |
                                            atlas::option::levels(srcField.levels()));

    // Target field initialization
    auto targetView = atlas::array::make_view<double, 2>(targetField);
    targetView.assign(0.0);

    // Interpolation from source field to target field
    interp_.execute(srcTmpField, targetField);

    // Redistribution from target field to destination field
    redistr_.execute(targetField, dstField);

    // Halo exchange
    dstField.haloExchange();
  }
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::Field & srcField,
                                        const atlas::Field & dstField) const {
  oops::Log::trace() << "AtlasInterpWrapper::executeAdjoint start "
                        "srcFieldName dstFieldName" <<
                        srcField.name() << " " << dstField.name() << std::endl;

  if (dstField.functionspace().type() == "PointCloud") {
    // PointCloud destination function space
    throw eckit::FunctionalityNotSupported("Adjoint not supported for PointCloud destination"
      " function space yet", Here());
  } else {
    // Other destination function space

    // Empty destination field setup
    atlas::Field dstTmp = dstField.functionspace().createField(dstField);

    // Copy of destination field
    const auto dstView = atlas::array::make_view<double, 2>(dstField);
    auto dstTmpView = atlas::array::make_view<double, 2>(dstTmp);
    for (atlas::idx_t t = 0; t < dstField.shape(0); ++t) {
      for (atlas::idx_t k = 0; k < dstField.shape(1); ++k) {
        dstTmpView(t, k) = dstView(t, k);
      }
    }

    // Halo exchange
    dstTmp.adjointHaloExchange();

    // Empty target field setup
    auto targetField = targetFunctionSpace_.createField<double>(
    atlas::option::name(dstField.name()) | atlas::option::levels(dstField.levels()));

    // Target field initialization
    auto targetView = atlas::array::make_view<double, 2>(targetField);
    targetView.assign(0.0);

    // Redistribution from destination field to target field
    inverseRedistr_.execute(dstTmp, targetField);

    // Source field initialization
    auto srcView = atlas::array::make_view<double, 2>(srcField);
    srcView.assign(0.0);

    // Adjoint interpolation target field to source field
    interp_.execute_adjoint(srcField, targetField);
  }

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
