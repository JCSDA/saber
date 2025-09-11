/*
 * (C) Copyright 2023-2024 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
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
  oops::Log::trace() << classname() << "::LinearVariableChange starting" << std::endl;

  // Local configuration
  LinearVariableChangeParameters params;
  params.deserialize(config);

  if (params.atlasFile.value() != boost::none) {
    // Check that input and output variables are present
    ASSERT(params.inputVariables.value() != boost::none);
    ASSERT(params.outputVariables.value() != boost::none);

    // Define output to input variables map
    ASSERT(params.inputVariables.value()->size() == params.outputVariables.value()->size());
    for (size_t jj = 0; jj < params.outputVariables.value()->size(); ++jj) {
      map_[(*params.outputVariables.value())[jj]] =
        (*params.inputVariables.value())[jj];
    }

    // Read multiplicative factor
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
    // TODO(AS): replace with setting up variables correctly
    util::readFieldSet(geom.getComm(),
                       geom.functionSpace(),
                       geom.variableSizes(*params.inputVariables.value()),
                       *params.inputVariables.value(),
                       conf,
                       multiplierFset_);
  }

  oops::Log::trace() << classname() << "::LinearVariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx,
                                       const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::changeVarTL starting" << std::endl;

  if (!multiplierFset_.empty()) {
    // Check that output variables are the same as in the map
    ASSERT(vars.size() == map_.size());
    for (auto const& it : map_) {
      ASSERT(vars.has(it.first));
    }

    // Get fieldset
    atlas::FieldSet fset;
    dx.toFieldSet(fset);

    // Multiply with the read fieldset
    util::multiplyFieldSets(fset, multiplierFset_);

    // Change variable name
    for (auto const& it : map_) {
      fset[it.second].rename(it.first);
    }

    // Back to increment
    dx.fromFieldSet(fset);
  }

  oops::Log::trace() << classname() << "::changeVarTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx,
                                              const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::changeVarInverseTL starting" << std::endl;

  if (!multiplierFset_.empty()) {
    // Check that output variables are the same as in the map
    ASSERT(vars.size() == map_.size());
    for (auto const& it : map_) {
      ASSERT(vars.has(it.second));
    }

    // Get fieldset
    atlas::FieldSet fset;
    dx.toFieldSet(fset);

    // Change variable name
    for (auto const& it : map_) {
      fset[it.first].rename(it.second);
    }

    // Divide by the read fieldset
    util::divideFieldSets(fset, multiplierFset_);

    // Back to increment
    dx.fromFieldSet(fset);
  }

  oops::Log::trace() << classname() << "::changeVarInverseTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx,
                                       const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::changeVarAD starting" << std::endl;

  if (!multiplierFset_.empty()) {
    // Check that output variables are the same as in the map
    ASSERT(vars.size() == map_.size());
    for (auto const& it : map_) {
      ASSERT(vars.has(it.second));
    }

    // Get fieldset
    atlas::FieldSet fset;
    dx.toFieldSet(fset);

    // Change variable name
    for (auto const& it : map_) {
      fset[it.first].rename(it.second);
    }

    // Multiply with the read fieldset
    util::multiplyFieldSets(fset, multiplierFset_);

    // Back to increment
    dx.fromFieldSet(fset);
  }

  oops::Log::trace() << classname() << "::changeVarAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                              const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::changeVarInverseAD starting" << std::endl;

  if (!multiplierFset_.empty()) {
    // Check that output variables are the same as in the map
    ASSERT(vars.size() == map_.size());
    for (auto const& it : map_) {
      ASSERT(vars.has(it.first));
    }

    // Get fieldset
    atlas::FieldSet fset;
    dx.toFieldSet(fset);

    // Divide by the read fieldset
    util::divideFieldSets(fset, multiplierFset_);

    // Change variable name
    for (auto const& it : map_) {
      fset[it.second].rename(it.first);
    }

    // Back to increment
    dx.fromFieldSet(fset);
  }

  oops::Log::trace() << classname() << "::changeVarInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
