/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/PressureParameters.h"

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

namespace saber {

oops::Variables GaussUVToGPParameters::mandatoryActiveVars() const {
  return oops::Variables({
      std::vector<std::string>{"eastward_wind",
                               "geostrophic_pressure_levels_minus_one",
                               "northward_wind"}});
}

oops::Variables GaussUVToGPParameters::activeInnerVars(const oops::Variables& outerVars) const {
  const int modelLevels = outerVars["eastward_wind"].getLevels();
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables vars;
  vars.push_back(oops::Variable{"eastward_wind", conf});
  vars.push_back(oops::Variable{"northward_wind", conf});
  return vars;
}

oops::Variables GpToHpParameters::mandatoryActiveVars() const {
  return oops::Variables({
      std::vector<std::string>{"geostrophic_pressure_levels_minus_one",
                               "hydrostatic_pressure_levels",
                               "unbalanced_pressure_levels_minus_one"}});
}

oops::Variables GpToHpParameters::activeInnerVars(const oops::Variables& outerVars) const {
  const int modelLevels = outerVars["hydrostatic_pressure_levels"].getLevels() - 1;
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables vars;
  vars.push_back(oops::Variable{"geostrophic_pressure_levels_minus_one", conf});
  vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
  return vars;
}

oops::Variables GpToHpm1Parameters::mandatoryActiveVars() const {
  return oops::Variables({
    std::vector<std::string>{"geostrophic_pressure_levels_minus_one",
                             "hydrostatic_pressure_levels_minus_one",
                             "unbalanced_pressure_levels_minus_one"}});
}

oops::Variables GpToHpm1Parameters::activeInnerVars(const oops::Variables& outerVars) const {
  const int modelLevels = outerVars["eastward_wind"].getLevels();
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables vars;
  vars.push_back(oops::Variable{"geostrophic_pressure_levels_minus_one", conf});
  vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
  return vars;
}

oops::Variables HydrostaticPressureParameters::mandatoryActiveVars() const {
  return oops::Variables({
    std::vector<std::string>{"eastward_wind",
                             "hydrostatic_pressure_levels",
                             "northward_wind",
                             "unbalanced_pressure_levels_minus_one"}});
}

oops::Variables HydrostaticPressureParameters::activeInnerVars(
    const oops::Variables& outerVars) const {
  const int modelLevels = outerVars["hydrostatic_pressure_levels"].getLevels() - 1;
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables vars;
  vars.push_back(oops::Variable{"eastward_wind", conf});
  vars.push_back(oops::Variable{"northward_wind", conf});
  vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
  return vars;
}

oops::Variables HydrostaticPressureParameters::intermediateTempVars(
    const oops::Variables& outerVars) const {
  if (outerVars.has("geostrophic_pressure_levels_minus_one")) {
    throw eckit::UserError("geostrophic_pressure_levels_minus_one is a "
        "temporary variable of mo_hydrostatic_pressure "
        " and should not be an outer variable of this block.",
        Here());
  }
  const int modelLevels = outerVars["hydrostatic_pressure_levels"].getLevels() - 1;
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables tempVars;
  tempVars.push_back(oops::Variable{"geostrophic_pressure_levels_minus_one", conf});
  return tempVars;
}

oops::Variables HydrostaticPressureMinusOneParameters::mandatoryActiveVars() const {
  return oops::Variables({
    std::vector<std::string>{"eastward_wind",
                             "hydrostatic_pressure_levels_minus_one",
                             "northward_wind",
                             "unbalanced_pressure_levels_minus_one"}});
}

oops::Variables HydrostaticPressureMinusOneParameters::activeInnerVars(
    const oops::Variables& outerVars) const {
  const int modelLevels = outerVars["eastward_wind"].getLevels();
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables vars;
  vars.push_back(oops::Variable{"eastward_wind", conf});
  vars.push_back(oops::Variable{"northward_wind", conf});
  vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
  return vars;
}

oops::Variables HydrostaticPressureMinusOneParameters::intermediateTempVars(
    const oops::Variables& outerVars) const {
  if (outerVars.has("geostrophic_pressure_levels_minus_one")) {
    throw eckit::UserError("geostrophic_pressure_levels_minus_one is a "
        "temporary variable of mo_hydrostatic_pressure "
        " and should not be an outer variable of this block.",
        Here());
  }
  const int modelLevels = outerVars["hydrostatic_pressure_levels_minus_one"].getLevels();
  eckit::LocalConfiguration conf;
  conf.set("levels", modelLevels);
  oops::Variables tempVars;
  tempVars.push_back(oops::Variable{"geostrophic_pressure_levels_minus_one", conf});
  return tempVars;
}

}  // namespace saber
