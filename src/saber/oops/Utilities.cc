/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "oops/base/FieldSet3D.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

oops::Variables getActiveVars(const SaberBlockParametersBase & params,
                              const oops::Variables & defaultVars) {
  oops::Log::trace() << "getActiveVars starting" << std::endl;
  oops::Variables activeVars_nomd;
  if (params.mandatoryActiveVars().size() == 0) {
    // No mandatory active variables for this block
    activeVars_nomd = params.activeVars.value().get_value_or(defaultVars);
  } else {
    // Block with mandatory active variables
    activeVars_nomd = params.activeVars.value().get_value_or(params.mandatoryActiveVars());
    ASSERT(params.mandatoryActiveVars() <= activeVars_nomd);
  }
  // Copy the variables that exist in defaultVars from defaultVars (they have metadata
  // associated with them)
  oops::Variables activeVars;
  for (auto & var : activeVars_nomd) {
    if (defaultVars.has(var.name())) {
      activeVars.push_back(defaultVars[var.name()]);
    } else {
      activeVars.push_back(var);
    }
  }
  return activeVars;
}

// -----------------------------------------------------------------------------

// Return inner variables as outer variables + inner active variables
// Can be used to help define innerVars_ member in SABER outer blocks
oops::Variables getUnionOfInnerActiveAndOuterVars(const SaberBlockParametersBase & params,
                                                  const oops::Variables & outerVars) {
  oops::Variables innerVars(outerVars);
  innerVars += params.activeInnerVars(outerVars);
  return innerVars;
}

// -----------------------------------------------------------------------------

// Return inner variables that are not outer variables
oops::Variables getInnerOnlyVars(const SaberBlockParametersBase & params,
                                 const oops::Variables & outerVars) {
  oops::Variables innerOnlyVars(getUnionOfInnerActiveAndOuterVars(params, outerVars));
  innerOnlyVars -= outerVars;
  return innerOnlyVars;
}

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi) {
  oops::Log::trace() << "setMPI starting" << std::endl;

  if (conf.has("mpi pattern")) {
    std::string mpiPattern = conf.getString("mpi pattern");
    util::seekAndReplace(conf, mpiPattern, std::to_string(mpi));
  }

  oops::Log::trace() << "setMPI done" << std::endl;
}

// -----------------------------------------------------------------------------

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::Variables & vars) {
  for (const auto& var : vars) {
    if (fset.has(var.name())) {
      throw eckit::UserError("Variable " + var.name() + " is already allocated in FieldSet.",
                             Here());
    }
  }
}

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::Variables & varsToAllocate,
                           const oops::Variables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace) {
  oops::Log::trace() << "allocateMissingFields starting" << std::endl;
  for (const auto& var : varsToAllocate) {
    if (!fset.has(var.name())) {
      oops::Log::info() << "Info     : Allocating " << var.name() << std::endl;
      auto field = functionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(varsWithLevels[var.name()].getLevels()));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      field.set_dirty(false);
      fset.add(field);
    }
  }
  oops::Log::trace() << "allocateMissingFields done" << std::endl;
}

// -----------------------------------------------------------------------------


}  // namespace saber
