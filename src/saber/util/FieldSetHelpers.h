/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <tuple>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

namespace util {
/// Return Field values as a function of distance to specific points.
/// Mostly relevant if the field is invariant by rotation around the point,
/// e.g., for Dirac tests of isotropic operators.
std::tuple<std::vector<std::vector<double>>,
           std::vector<std::vector<double>>,
           std::vector<size_t>,
           std::vector<size_t>>
sortBySeparationDistance(const eckit::mpi::Comm &,
                         const atlas::FunctionSpace &,
                         const atlas::FieldSet &,
                         const atlas::FieldSet &,
                         const double,
                         const bool);
}  // namespace util
