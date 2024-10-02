/*
 * (C) Crown Copyright 2021-2024, Met Office
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/AtlasInterpWrapper.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "eckit/exception/Exceptions.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

AtlasInterpWrapper::AtlasInterpWrapper(const atlas::grid::Partitioner & srcPartitioner,
                                       const atlas::FunctionSpace & srcFspace,
                                       const atlas::Grid & dstGrid,
                                       const atlas::FunctionSpace & dstFspace,
                                       const std::string & interpType)
  : targetFspace_(), interp_(), redistr_(), inverseRedistr_() {
  oops::Log::trace() << classname() << "::AtlasInterpWrapper starting" << std::endl;

  // Get or compute source mesh
  atlas::Mesh srcMesh;
  if (srcFspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns fs(srcFspace);
    srcMesh = atlas::MeshGenerator("structured").generate(fs.grid(), srcPartitioner);
  } else if (srcFspace.type() == "NodeColumns") {
    const atlas::functionspace::NodeColumns fs(srcFspace);
    srcMesh = fs.mesh();
  } else {
    throw eckit::FunctionalityNotSupported(srcFspace.type()
      + " source function space not supported yet", Here());
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
    throw eckit::FunctionalityNotSupported(dstFspace.type()
      + " destination function space type not supported yet", Here());
    targetFspace_ = atlas::FunctionSpace();
  }

  // Interpolation
  atlas::util::Config interpConfig;
  if (srcFspace.type() == "StructuredColumns") {
    // StructuredColumns
    interpConfig.set("type",
                     interpType == "" ? "structured-bilinear" : interpType);
  } else if (srcFspace.type() == "NodeColumns") {
    // NodeColumns
    interpConfig.set("type",
                     interpType == "" ? "unstructured-bilinear-lonlat" : interpType);
  } else {
    throw eckit::FunctionalityNotSupported(srcFspace.type()
      + " source function space type not supported yet", Here());
  }
  interpConfig.set("adjoint", "true");
  interp_ = atlas::Interpolation(interpConfig, srcFspace, targetFspace_);

  // Redistribution
  redistr_ = atlas::Redistribution(targetFspace_, dstFspace);

  // Inverse redistribution
  inverseRedistr_ = atlas::Redistribution(dstFspace, targetFspace_);


  oops::Log::trace() << classname() << "::AtlasInterpWrapper done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::Field & srcField,
                                 atlas::Field & dstField) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  // Empty source field setup
  atlas::Field srcTmpField = srcField.functionspace().createField(srcField);

  // Copy of source field (this includes halo points; these were updated before calling execute)
  const auto srcView = atlas::array::make_view<double, 2>(srcField);
  auto srcTmpView = atlas::array::make_view<double, 2>(srcTmpField);
  for (atlas::idx_t t = 0; t < srcField.shape(0); ++t) {
    for (atlas::idx_t k = 0; k < srcField.shape(1); ++k) {
      srcTmpView(t, k) = srcView(t, k);
    }
  }

  // Empty target field setup
  auto targetField = targetFspace_.createField<double>(
    atlas::option::name(dstField.name()) | atlas::option::levels(srcField.shape(1)));

  // Target field initialization
  auto targetView = atlas::array::make_view<double, 2>(targetField);
  targetView.assign(0.0);

  // Interpolation from source field to target field
  interp_.execute(srcTmpField, targetField);

  // Redistribution from target field to destination field
  redistr_.execute(targetField, dstField);

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::Field & srcField,
                                        const atlas::Field & dstField) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

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

  // Empty target field setup
  auto targetField = targetFspace_.createField<double>(
  atlas::option::name(dstField.name()) | atlas::option::levels(dstField.shape(1)));

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

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::execute(const atlas::FieldSet & srcFieldSet,
                                 atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  srcFieldSet.haloExchange();
  for (auto & srcField : srcFieldSet) {
    execute(srcField, targetFieldSet[srcField.name()]);
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void AtlasInterpWrapper::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                        const atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  for (auto & srcField : srcFieldSet) {
    executeAdjoint(srcField, targetFieldSet[srcField.name()]);
  }
  srcFieldSet.adjointHaloExchange();
  srcFieldSet.set_dirty();

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace saber
