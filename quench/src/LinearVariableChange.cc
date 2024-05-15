/*
 * (C) Copyright 2023  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/LinearVariableChange.h"

#include <ostream>
#include <string>

#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

#include "src/LinearVariableChangeParameters.h"

namespace quench {

// -----------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const eckit::Configuration & config) {
  // Local configuration
  LinearVariableChangeParameters params;
  params.deserialize(config);

  if (params.variables.value() != boost::none) {
    // Read multiplicative factor
    ASSERT(params.atlasFile.value() != boost::none);
    eckit::LocalConfiguration conf(*params.atlasFile.value());

    // Get number of MPI tasks and OpenMP threads
    std::string mpi(std::to_string(geom.getComm().size()));
    std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
    {
      omp = std::to_string(omp_get_num_threads());
    }
#endif

    // Replace patterns
    util::seekAndReplace(conf, "_MPI_", mpi);
    util::seekAndReplace(conf, "_OMP_", omp);

    // Read fieldset
    util::readFieldSet(geom.getComm(),
                       geom.functionSpace(),
                       geom.variableSizes(*params.variables.value()),
                       params.variables.value()->variables(),
                       conf,
                       fset_);
  }
}

// -----------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars) const {
  if (!fset_.empty()) {
    atlas::FieldSet fset;
    dx.toFieldSet(fset);
    util::multiplyFieldSets(fset, fset_);
    dx.fromFieldSet(fset);
  }
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx, const oops::Variables & vars) const {
  if (!fset_.empty()) {
    atlas::FieldSet fset;
    dx.toFieldSet(fset);
    util::divideFieldSets(fset, fset_);
    dx.fromFieldSet(fset);
  }
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars) const {
  if (!fset_.empty()) {
    atlas::FieldSet fset;
    dx.toFieldSet(fset);
    util::multiplyFieldSets(fset, fset_);
    dx.fromFieldSet(fset);
  }
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx, const oops::Variables & vars) const {
  if (!fset_.empty()) {
    atlas::FieldSet fset;
    dx.toFieldSet(fset);
    util::divideFieldSets(fset, fset_);
    dx.fromFieldSet(fset);
  }
}

// -----------------------------------------------------------------------------

}  // namespace quench
