/*
 * (C) Crown Copyright 2021-2024, Met Office
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/interpolation/VectorFieldMetadata.h"
#include "saber/util/defines.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/for_each.h"

#include "eckit/exception/Exceptions.h"
// -----------------------------------------------------------------------------

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

AtlasInterpWrapper::AtlasInterpWrapper(const atlas::grid::Partitioner & srcPartitioner,
                                       const atlas::FunctionSpace & srcFspace,
                                       const atlas::Grid & dstGrid,
                                       const atlas::FunctionSpace & dstFspace,
                                       const std::string & interpType,
                                       const bool includingVectorInterpolation)
  : targetFspace_(), interp_(), redistr_(), inverseRedistr_(),
    includingVectorInterpolation_(includingVectorInterpolation) {
  oops::Log::trace() << classname() << "::AtlasInterpWrapper starting" << std::endl;

  const auto startMPIComm = eckit::mpi::comm().name();

  // Get or compute source mesh
  atlas::Mesh srcMesh;
  if (srcFspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns fs(srcFspace);

    if constexpr (!SABER_ATLAS_VERSION_46_OR_GREATER) {
      eckit::mpi::setCommDefault(startMPIComm);
    }

    srcMesh = atlas::MeshGenerator("structured").generate(fs.grid(), srcPartitioner);
  } else if (srcFspace.type() == "NodeColumns") {
    const atlas::functionspace::NodeColumns fs(srcFspace);
    srcMesh = fs.mesh();
  } else {
    throw eckit::NotImplemented(srcFspace.type()
      + " source function space not supported yet", Here());
  }

  if constexpr (!SABER_ATLAS_VERSION_46_OR_GREATER) {
    eckit::mpi::setCommDefault(startMPIComm);
  }

  // Get target partitioner from source mesh
  const atlas::grid::MatchingMeshPartitioner targetPartitioner(srcMesh);

  // Create target function space
  if (dstFspace.type() == "StructuredColumns") {
    // StructuredColumns
    const atlas::grid::Distribution distribution(dstGrid, targetPartitioner);
    targetFspace_ = atlas::functionspace::StructuredColumns(dstGrid, distribution);
  } else if (dstFspace.type() == "NodeColumns") {
    // NodeColumns
    if (dstGrid.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      const atlas::Mesh targetMesh = atlas::MeshGenerator("cubedsphere_dual").generate(dstGrid,
        targetPartitioner);
      targetFspace_ = atlas::functionspace::CubedSphereNodeColumns(targetMesh);
    } else {
      // Other NodeColumns (TODO: not working!)
      const atlas::Mesh targetMesh = atlas::MeshGenerator("delaunay").generate(dstGrid,
        targetPartitioner);
      targetFspace_ = atlas::functionspace::NodeColumns(targetMesh);
    }
  } else {
    throw eckit::NotImplemented(dstFspace.type()
      + " destination function space type not supported yet", Here());
  }

  // Ensure comm definitely set correctly prior to interpolation.
  if constexpr (!SABER_ATLAS_VERSION_46_OR_GREATER) {
    eckit::mpi::setCommDefault(startMPIComm);
  }

  // Interpolation
  atlas::util::Config interpScalarConfig;
  if (srcFspace.type() == "StructuredColumns") {
    // StructuredColumns
    interpScalarConfig.set("type",
                     interpType == "" ? "structured-bilinear" : interpType);
    interpScalarConfig.set("halo",
                     interpType.find("cubic") != interpType.npos ? 2 : 1);
  } else if (srcFspace.type() == "NodeColumns") {
    // NodeColumns
    interpScalarConfig.set("type",
                     interpType == "" ? "unstructured-bilinear-lonlat" : interpType);
    interpScalarConfig.set("halo", 1);
  } else {
    throw eckit::NotImplemented(srcFspace.type()
      + " source function space type not supported yet", Here());
  }

  interpScalarConfig.set("adjoint", "true");

  atlas::util::Config interpConfig;
  if (includingVectorInterpolation_) {
    static const auto sphericalVector =
      atlas::option::type("spherical-vector") | atlas::util::Config("adjoint", true);
    interpConfig = sphericalVector | atlas::util::Config("scheme", interpScalarConfig);

  } else {
    interpConfig = interpScalarConfig;
  }

  oops::Log::info() << "interpConfig is: " << interpConfig << std::endl;

  interp_ = atlas::Interpolation(interpConfig, srcFspace, targetFspace_);

  // Redistribution
  redistr_ = atlas::Redistribution(targetFspace_, dstFspace);

  // Inverse redistribution
  inverseRedistr_ = atlas::Redistribution(dstFspace, targetFspace_);

  oops::Log::trace() << classname() << "::AtlasInterpWrapper done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::FieldSet & srcFieldSet,
                                 atlas::FieldSet & dstFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  srcFieldSet.haloExchange();

  atlas::FieldSet tmpSrcFieldSet = atlas::FieldSet{};
  std::vector<size_t> variableSizes;

  for (auto& srcField : srcFieldSet) {
    atlas::Field srcTmpField =
      srcField.functionspace().createField<double>(
                                 atlas::option::name(srcField.name()) |
                                 atlas::option::levels(srcField.shape(1)));

    util::for_each_value([](const double rhs, double & lhs) { lhs = rhs; },
                         srcField,
                         srcTmpField);
    srcTmpField.set_dirty();

    tmpSrcFieldSet.add(srcTmpField);

    variableSizes.push_back(srcField.levels());
  }

  atlas::FieldSet targetFieldSet = util::createFieldSet(targetFspace_, variableSizes,
                                                        dstFieldSet.field_names(), 0.0);

  if (includingVectorInterpolation_) {
    appendVectorFieldMeta(tmpSrcFieldSet);
    appendVectorFieldMeta(targetFieldSet);
  }

  interp_.execute(tmpSrcFieldSet, targetFieldSet);

  redistr_.execute(targetFieldSet, dstFieldSet);

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                        const atlas::FieldSet & dstFieldSet) const {
    oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

    atlas::FieldSet tmpDstFieldSet = atlas::FieldSet{};
    std::vector<size_t> variableSizes;

    for (auto& dstField : dstFieldSet) {
      atlas::Field tmpDstField =
        dstField.functionspace().createField<double>(
                                   atlas::option::name(dstField.name()) |
                                   atlas::option::levels(dstField.shape(1)));

      util::for_each_value([](const double rhs, double & lhs) { lhs = rhs; },
                           dstField,
                           tmpDstField);
      tmpDstField.set_dirty();

      tmpDstFieldSet.add(tmpDstField);

      variableSizes.push_back(dstField.levels());

      auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[dstField.name()]);
      srcView.assign(0.0);
    }

    atlas::FieldSet targetFieldSet = util::createFieldSet(targetFspace_, variableSizes,
                                                          dstFieldSet.field_names(), 0.0);

    if (includingVectorInterpolation_) {
      appendVectorFieldMeta(srcFieldSet);
      appendVectorFieldMeta(targetFieldSet);
    }

    inverseRedistr_.execute(tmpDstFieldSet, targetFieldSet);

    interp_.execute_adjoint(srcFieldSet, targetFieldSet);

    srcFieldSet.adjointHaloExchange();
    srcFieldSet.set_dirty();

    oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace saber
