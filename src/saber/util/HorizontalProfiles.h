/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

namespace util {
/// Return Field values as a function of distance to specific points.
/// Mostly relevant if the field is invariant by rotation around the point,
/// e.g., for Dirac tests of isotropic operators.
/// Note the sorting by distance includes a tolerance to obtain portable
/// behavior -- however this means the sorted points may not be in strictly
/// monotonically-increasing distance (at the 1e-13 level).
std::tuple<std::vector<std::vector<double>>,
           std::vector<std::vector<double>>,
           std::vector<double>,
           std::vector<double>,
           std::vector<size_t>,
           std::vector<size_t>>
sortBySeparationDistance(const eckit::mpi::Comm & comm,
                         const atlas::FunctionSpace & fspace,
                         const atlas::FieldSet & dataFset,
                         const atlas::FieldSet & diracFset,
                         const double maxLength,
                         const bool removeDuplicates);

// Write horizontal covariances to file
void  write_1d_covariances(const eckit::mpi::Comm & comm,
                           const std::vector<std::vector<double>> & distances,
                           const std::vector<std::vector<double>> & covariances,
                           const std::vector<double> & lons,
                           const std::vector<double> & lats,
                           const std::vector<size_t> & levs,
                           const std::vector<size_t> & fieldIndexes,
                           const std::vector<std::string> & names,
                           const std::string & filePath);
}  // namespace util
