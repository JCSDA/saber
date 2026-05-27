/*
 * (C) Copyright 2026- Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

namespace util {

// -----------------------------------------------------------------------------

void randomCtlVec(const eckit::mpi::Comm &,
                  const std::vector<int> &,
                  std::vector<double> &);

// -----------------------------------------------------------------------------

void randomCtlVec(const eckit::mpi::Comm &,
                  const std::vector<int> &,
                  atlas::Field &);

// -----------------------------------------------------------------------------

}  // namespace util

