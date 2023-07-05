/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/Utilities.h"

// -----------------------------------------------------------------------------
namespace saber {

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
    return oops::Variables(varsconf, varsStrings);
  } else {
    return activeVars;
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
