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
  oops::Variables activeVars;
  if (params.mandatoryActiveVars().size() == 0) {
    // No mandatory active variables for this block
    activeVars = params.activeVars.value().get_value_or(defaultVars);
  } else {
    // Block with mandatory active variables
    activeVars = params.activeVars.value().get_value_or(params.mandatoryActiveVars());
    ASSERT(params.mandatoryActiveVars() <= activeVars);
  }
  if (activeVars.variablesMetaData().empty()) {
    atlas::util::Config defvarsconf(defaultVars.variablesMetaData());
    atlas::util::Config varsconf;
    std::vector<std::string> varsStrings(activeVars.variables());
    for (const std::string & var : activeVars.variables()) {
      varsconf = varsconf | atlas::util::Config(var, defvarsconf.getSubConfiguration(var));
    }
    activeVars = oops::Variables(varsconf, varsStrings);
  }
  return activeVars;
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

void allocateFields(oops::FieldSet3D & fset,
                    const oops::Variables & varsToAllocate,
                    const oops::Variables & varsWithLevels,
                    const atlas::FunctionSpace & functionSpace,
                    const bool haloExchange) {
  oops::Log::trace() << "allocateFields starting" << std::endl;
  for (const auto& var : varsToAllocate.variables()) {
    if (!fset.has(var)) {
      oops::Log::test() << "Allocating " << var << std::endl;
      auto field = functionSpace.createField<double>(
                atlas::option::name(var) |
                atlas::option::levels(varsWithLevels.getLevels(var)));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      if (haloExchange) {
        field.haloExchange();
      }
      fset.add(field);
    }
  }
  oops::Log::trace() << "allocateFields done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
