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

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class AirTemperatureParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperatureParameters, SaberBlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "air_temperature",
    "exner_levels_minus_one",
    "potential_temperature"});}
};

// -----------------------------------------------------------------------------

class MoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistIncrOpParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<std::string> mio_file{"moisture incrementing operator file", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "air_temperature",
    "mass_content_of_cloud_ice_in_atmosphere_layer",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "qt",
    "specific_humidity"});}
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
};


}  // namespace saber
