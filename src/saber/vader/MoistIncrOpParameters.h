/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class AirTemperatureParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperatureParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({
       "air_temperature",
       "exner_levels_minus_one",
       "potential_temperature"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"potential_temperature",
                          "exner_levels_minus_one"});
    const int modelLevels = outerVars.getLevels("air_temperature");
    vars.addMetaData("potential_temperature", "levels", modelLevels);
    vars.addMetaData("exner_levels_minus_one", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"air_temperature"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

class MoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistIncrOpParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<std::string> mio_file{"moisture incrementing operator file", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "air_temperature",
    "mass_content_of_cloud_ice_in_atmosphere_layer",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "qt",
    "specific_humidity"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"air_temperature",
                          "qt"});
    const int modelLevels = outerVars.getLevels("specific_humidity");
    vars.addMetaData("air_temperature", "levels", modelLevels);
    vars.addMetaData("qt", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"mass_content_of_cloud_ice_in_atmosphere_layer",
                          "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
                          "specific_humidity"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

class SuperMoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperMoistIncrOpParameters, SaberBlockParametersBase)

 public:
  AirTemperatureParameters airTemperature{this};
  MoistIncrOpParameters moistIncrOp{this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "exner_levels_minus_one",
    "potential_temperature",
    "mass_content_of_cloud_ice_in_atmosphere_layer",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "qt",
    "specific_humidity"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"exner_levels_minus_one",
                          "potential_temperature",
                          "qt"});
    const int modelLevels = outerVars.getLevels("specific_humidity");
    vars.addMetaData("exner_levels_minus_one", "levels", modelLevels);
    vars.addMetaData("potential_temperature", "levels", modelLevels);
    vars.addMetaData("qt", "levels", modelLevels);
    return vars;
  }

  // activeOuterVars() not needed in this super block.
  // It would have contained "mass_content_of_cloud_ice_in_atmosphere_layer",
  // "mass_content_of_cloud_liquid_water_in_atmosphere_layer" and "specific_humidity".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
    oops::Variables tempVars({"air_temperature"});
    if (outerVars.has("air_temperature")) {
      throw eckit::UserError("air_temperature is a temporary variable of mo_super_mio"
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels = outerVars.getLevels("specific_humidity");
    tempVars.addMetaData("air_temperature", "levels", modelLevels);
    return tempVars;
  }
};


}  // namespace saber
