/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

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
    return oops::Variables({std::vector<std::string>{
       "air_temperature",
       "dimensionless_exner_function_levels_minus_one",
       "air_potential_temperature"}});
  }

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({
       "height_above_mean_sea_level",
       "height_above_mean_sea_level_levels",
       "dimensionless_exner_function_levels_minus_one",
       "air_potential_temperature"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["air_temperature"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"air_potential_temperature", conf});
    vars.push_back({"dimensionless_exner_function_levels_minus_one", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["air_temperature"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class MoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistIncrOpParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<std::string> mio_file{"moisture incrementing operator file", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "air_temperature",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "qt",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"}});}

  const oops::Variables mandatoryStateVars() const override {return oops::Variables({
    "liquid_cloud_volume_fraction_in_atmosphere_layer",
    "ice_cloud_volume_fraction_in_atmosphere_layer",
    "qsat", "dlsvpdT", "rht"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels =
              outerVars["water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"air_temperature", conf});
    vars.push_back({"qt", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars(
      {outerVars["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"],
       outerVars["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"],
       outerVars["water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"]});
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
    std::vector<std::string>{
    "dimensionless_exner_function_levels_minus_one",
    "air_potential_temperature",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "qt",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"}});}

  // combined variables for AirTemperature and MoistIncrOp
  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({
       "height_above_mean_sea_level",
       "height_above_mean_sea_level_levels",
       "dimensionless_exner_function_levels_minus_one",
       "air_potential_temperature",
       "liquid_cloud_volume_fraction_in_atmosphere_layer",
       "ice_cloud_volume_fraction_in_atmosphere_layer",
       "qsat", "dlsvpdT", "rht"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels =
              outerVars["water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"dimensionless_exner_function_levels_minus_one", conf});
    vars.push_back({"air_potential_temperature", conf});
    vars.push_back({"qt", conf});
    return vars;
  }

  // activeOuterVars() not needed in this super block.
  // It would have contained "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
  // "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water" and
  // "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water" -> AKA "specific_humidity".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
    if (outerVars.has("air_temperature")) {
      throw eckit::UserError("air_temperature is a temporary variable of mo_super_mio"
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels =
              outerVars["water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables tempVars;
    tempVars.push_back({"air_temperature", conf});
    return tempVars;
  }
};

}  // namespace saber
