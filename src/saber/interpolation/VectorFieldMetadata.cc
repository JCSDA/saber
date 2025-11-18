/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string_view>
#include <unordered_map>

#include "atlas/library/version.h"
#include "atlas/option/Options.h"
#include "atlas/util/Config.h"

#include "saber/interpolation/VectorFieldMetadata.h"


namespace saber {
namespace interpolation {

namespace {

// Retrieves the vector field configuration for a given variable name, if there is one.
atlas::util::Config getVectorFieldOpt(const std::string_view varName) {
#ifdef ATLAS_SUPPORTS_SPHERICAL_VECTOR_INTERP_META
  constexpr size_t xComponent = 0;
  constexpr size_t yComponent = 1;

  if (varName == "eastward_wind") {
    return atlas::option::vector_component("wind", xComponent);
  } else if (varName == "northward_wind") {
    return atlas::option::vector_component("wind", yComponent);
  } else if (varName == "eastward_wind_at_10m") {
    return atlas::option::vector_component("wind_at_10m", xComponent);
  } else if (varName == "northward_wind_at_10m") {
    return atlas::option::vector_component("wind_at_10m", yComponent);
  }

#endif  // ATLAS_SUPPORTS_SPHERICAL_VECTOR_INTERP_META

  // Empty config - don't add any vector config to this field.
  return atlas::util::Config();
}

}  // namespace

// ------------------------------------------------------------------------------------------------

void appendVectorFieldMeta(atlas::FieldSet& fset) {
  // If supported, use spherical vector implementation.
#ifdef ATLAS_SUPPORTS_SPHERICAL_VECTOR_INTERP_META
  for (auto& field : fset) {
    auto opt = getVectorFieldOpt(field.name());
    if (!opt.empty()) {   // If this field is a component of a vector field.
      field.metadata().set(opt);
    }
  }
#else
  // Otherwise, 'old-school' method.
  if (fset.has("eastward_wind") && fset.has("northward_wind")) {
    fset["eastward_wind"].metadata().set("vector_field_name", "wind");
    fset["northward_wind"].metadata().set("vector_field_name", "wind");
  }
#endif  // ATLAS_SUPPORTS_SPHERICAL_VECTOR_INTERP_META
}

}  // namespace interpolation
}  // namespace saber
