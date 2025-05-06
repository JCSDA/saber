/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

namespace atlas {
class FunctionSpace;
}  // namespace atlas

namespace saber::gsi {

// Extract from an atlas FunctionSpace the data that will enable comparing it to the
// GSI grid in Fortran code.
//
// If the FunctionSpace is StructuredColumns with a regular checkerboard partitioning
// over MPI tasks, then the grid's nx/ny and its x/y coordinates will be extracted and
// placed into the return vector.
//
// If not, the function will throw an exception as the data to compare with GSI can't
// even be extracted.
std::vector<double> functionspaceToGridChecks(const atlas::FunctionSpace & fs);

}  // namespace saber::gsi
