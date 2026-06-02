/*
 * (C) Crown Copyright 2026-, Met Office
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "atlas/grid/Distribution.h"
#include "atlas/meshgenerator/MeshGenerator.h"

#include "eckit/exception/Exceptions.h"

#include "saber/interpolation/SMVInterpWrapper.h"
#include "saber/interpolation/VectorFieldMetadata.h"
#include "saber/util/defines.h"

#include "oops/util/Logger.h"

namespace saber {
namespace interpolation {

SMVInterpWrapper::SMVInterpWrapper(
    const atlas::functionspace::NodeColumns& CSFunctionSpace,
    const atlas::functionspace::StructuredColumns& gaussFunctionSpace,
    const atlas::Grid& gaussGrid, const bool includingVectorInterpolation,
    const bool normalizingInterpolation)
    : gaussFunctionSpace_(gaussFunctionSpace),
      includingVectorInterpolation_(includingVectorInterpolation),
      normalizingInterpolation_(normalizingInterpolation) {
  if constexpr (!SABER_ATLAS_VERSION_46_OR_GREATER) {
    const std::string errMsg =
        "Atlas 0.46.0 or newer is required for interpolating using the SMV "
        "interpolator. Aborting.";
    eckit::Exception(errMsg, Here());
  } else {
    oops::Log::trace() << "::SMVInterpWrapper\n";
  }

  // Construct intermediate node columns function space with Gauss mesh
  const atlas::grid::MatchingMeshPartitioner CSPartitioner(
      CSFunctionSpace.mesh(), atlas::option::type("cubedsphere"));
  const atlas::grid::Distribution gaussToCSDistribution(gaussGrid, CSPartitioner);
  const atlas::Mesh structuredGaussMesh =
      atlas::StructuredMeshGenerator()
          .generate(gaussGrid, gaussToCSDistribution);

  gaussToCSFunctionSpace_ =
      atlas::functionspace::NodeColumns(structuredGaussMesh);

  // Set up interpolation object
  atlas::util::Config interpScalarConfig;
  interpScalarConfig.set("type", "spherical-mean-value");
  interpScalarConfig.set("normalised", normalizingInterpolation);
  interpScalarConfig.set("adjoint", false);

  atlas::util::Config interpConfig;
  if (includingVectorInterpolation_) {
    interpConfig = atlas::option::type("spherical-vector") |
                   atlas::util::Config("scheme", interpScalarConfig);
  } else {
    interpConfig = interpScalarConfig;
  }

  interpolation_ =
      atlas::Interpolation(interpConfig, CSFunctionSpace, gaussToCSFunctionSpace_);
  redistribution_ =
      atlas::Redistribution(gaussToCSFunctionSpace_, gaussFunctionSpace);
}

void SMVInterpWrapper::doInterpolation(
    const std::vector<atlas::util::Config>& fieldConfigs,
    const atlas::FieldSet& srcFieldSet,
    atlas::FieldSet& newFieldSet) const {
  oops::Log::trace() << "::SMVInterpWrapper::doInterpolation\n";
  srcFieldSet->haloExchange();

  // Set up fields for CS and gaussToCS function spaces
  atlas::FieldSet gaussToCSFieldSet;
  atlas::FieldSet tmpSrcFieldSet;

  for (const atlas::util::Config& fieldConfig : fieldConfigs) {
    const std::string fieldName = fieldConfig.getString("name");

    atlas::Field gaussToCSField =
        gaussToCSFunctionSpace_->createField<double>(fieldConfig);
    atlas::Field tmpSrcField =
        srcFieldSet[fieldName].functionspace().createField<double>(fieldConfig);

    atlas::array::make_view<double, 2>(gaussToCSField).assign(0.0);
    const auto srcView =
        atlas::array::make_view<double, 2>(srcFieldSet[fieldName]);
    auto srcTmpView = atlas::array::make_view<double, 2>(tmpSrcField);
    for (atlas::idx_t t = 0; t < srcFieldSet[fieldName].shape(0); ++t) {
      for (atlas::idx_t k = 0; k < srcFieldSet[fieldName].shape(1); ++k) {
        srcTmpView(t, k) = srcView(t, k);
      }
    }

    gaussToCSFieldSet.add(gaussToCSField);
    tmpSrcFieldSet.add(tmpSrcField);
  }

  if (includingVectorInterpolation_) {
    appendVectorFieldMeta(gaussToCSFieldSet);
    appendVectorFieldMeta(tmpSrcFieldSet);
  }

  // Do interpolation and redistribution
  interpolation_->execute(tmpSrcFieldSet, gaussToCSFieldSet);

  for (const atlas::util::Config& fieldConfig : fieldConfigs) {
    atlas::Field gaussField = gaussFunctionSpace_->createField<double>(fieldConfig);
    newFieldSet.add(gaussField);
  }

  redistribution_->execute(gaussToCSFieldSet, newFieldSet);
}

}  // namespace interpolation
}  // namespace saber
