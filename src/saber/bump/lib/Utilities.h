/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"

namespace bump_lib {

// -----------------------------------------------------------------------------
// ConfigFunctions
// -----------------------------------------------------------------------------

bool isVector(const eckit::Configuration &);

bool isSubConfig(const eckit::Configuration &);

bool isFinal(const eckit::Configuration &);

void seekAndReplace(eckit::LocalConfiguration &,
                    const std::string &,
                    const std::string &);
void seekAndReplace(eckit::LocalConfiguration &,
                    const std::string &,
                    const size_t &,
                    const size_t &);

eckit::LocalConfiguration mergeConfigs(const eckit::Configuration &,
                                       const eckit::Configuration &);

// -----------------------------------------------------------------------------
// FieldSetHelpers
// -----------------------------------------------------------------------------

void shareFields(const atlas::FieldSet &, atlas::FieldSet &);

std::string getGridUid(const atlas::FunctionSpace &);

std::string getGridUid(const atlas::FieldSet &);

void readFieldSet(const eckit::mpi::Comm &,
                  const atlas::FunctionSpace &,
                  const std::vector<size_t> &,
                  const std::vector<std::string> &,
                  const eckit::Configuration &,
                  atlas::FieldSet &);

void writeFieldSet(const eckit::mpi::Comm &,
                   const eckit::Configuration &,
                   const atlas::FieldSet &);

// -----------------------------------------------------------------------------
// FieldSetOperations
// -----------------------------------------------------------------------------

double normFieldSet(const atlas::FieldSet &,
                    const std::vector<std::string> &,
                    const eckit::mpi::Comm &);

// -----------------------------------------------------------------------------

}  // namespace bump_lib
