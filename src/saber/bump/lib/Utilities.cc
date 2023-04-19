/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/lib/Utilities.h"

#include <omp.h>

#include <string>

#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

namespace bump_lib {

// -----------------------------------------------------------------------------
// ConfigFunctions from OOPS
// -----------------------------------------------------------------------------

bool isVector(const eckit::Configuration & config) {return util::isVector(config);}

bool isSubConfig(const eckit::Configuration & config) {return util::isSubConfig(config);}

bool isFinal(const eckit::Configuration & config) {return util::isFinal(config);}

void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                    const std::string & value)
  {util::seekAndReplace(config, pattern, value);}
void seekAndReplace(eckit::LocalConfiguration & config, const std::string & pattern,
                    const size_t & count, const size_t & zpad)
  {util::seekAndReplace(config, pattern, count, zpad);}

eckit::LocalConfiguration mergeConfigs(const eckit::Configuration & config1,
                                       const eckit::Configuration & config2)
  {return util::mergeConfigs(config1, config2);}

// -----------------------------------------------------------------------------
// FieldSetHelpers from OOPS
// -----------------------------------------------------------------------------

void shareFields(const atlas::FieldSet & otherFset, atlas::FieldSet & fset)
  {util::shareFields(otherFset, fset);}

std::string getGridUid(const atlas::FunctionSpace & fspace) {return util::getGridUid(fspace);}

std::string getGridUid(const atlas::FieldSet & fset) {return util::getGridUid(fset);}

void readFieldSet(const eckit::mpi::Comm & comm,
                  const atlas::FunctionSpace & fspace,
                  const std::vector<size_t> & variableSizes,
                  const std::vector<std::string> & vars,
                  const eckit::Configuration & config,
                  atlas::FieldSet & fset)
  {return util::readFieldSet(comm, fspace, variableSizes, vars, config, fset);}

void writeFieldSet(const eckit::mpi::Comm & comm,
                   const eckit::Configuration & config,
                   const atlas::FieldSet & fset)
  {return util::writeFieldSet(comm, config, fset);}

// -----------------------------------------------------------------------------
// FieldSetOperations from OOPS
// -----------------------------------------------------------------------------

double normFieldSet(const atlas::FieldSet & fset,
                    const std::vector<std::string> & vars,
                    const eckit::mpi::Comm & comm)
  {return util::normFieldSet(fset, vars, comm);}

// -----------------------------------------------------------------------------

}  // namespace bump_lib
