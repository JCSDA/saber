/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include "atlas/field/FieldSet.h"

namespace saber {
namespace interpolation {

// Append metadata to fields for vector field interpolation.
void appendVectorFieldMeta(atlas::FieldSet& fset);

}  // namespace interpolation
}  // namespace saber
