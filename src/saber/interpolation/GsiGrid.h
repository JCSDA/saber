/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

namespace atlas {
  class Grid;
  class FunctionSpace;
  class FieldSet;
  class RegularGrid;
}

namespace eckit {
  class Configuration;
  namespace mpi {
    class Comm;
  }
}

namespace saber {
namespace interpolation {

/// @brief  Config key selecting the GSI-matching grid
extern const std::string GsiGridKey;
/// @brief  Config key selecting the GSI-matching MPI partitioner
///         (must be used with GsiGridKey)
extern const std::string GsiPartitionerKey;

/// @brief  Compute a south-to-north "checkerboard" MPI partition of a regular grid,
///         matching GSI's layout (GSI orders points/partitions south-to-north,
///         atlas north-to-south). `nbands` must divide `ntasks`.
std::vector<int> computeS2NCheckerboardPartition(const atlas::RegularGrid & rg,
                                                 const int ntasks, const int nbands);

/// @brief  Build an atlas grid + function space replicating the GSI grid and MPI
///         partition (but not GSI's per-task grid-point ordering).
void setupGsiMatchingGrid(const eckit::Configuration & config,
                          const eckit::mpi::Comm & comm,
                          atlas::Grid & grid,
                          atlas::FunctionSpace & functionSpace,
                          atlas::FieldSet & fieldSet);

}  // namespace interpolation
}  // namespace saber
