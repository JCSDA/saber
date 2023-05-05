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

class GaussUVToGPParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPParameters, SaberBlockParametersBase)
 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "eastward_wind",
    "geostrophic_pressure_levels_minus_one",
    "northward_wind"});}
};

class GpToHpCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GpToHpCovarianceParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::RequiredParameter<int> covariance_nlat{"number of covariance latitude rings", this};
  oops::Parameter<int> gp_regression_bins{"gp regression bins", "gP regression bins", 18, this};
};

// -----------------------------------------------------------------------------

class GpToHpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GpToHpParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<GpToHpCovarianceParameters>
    gptohpcovarianceparams{"covariance data", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "geostrophic_pressure_levels_minus_one",
    "hydrostatic_pressure_levels",
    "unbalanced_pressure_levels_minus_one"});}
};


class HydrostaticPressureParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticPressureParameters, SaberBlockParametersBase)
 public:
  GaussUVToGPParameters gaussUVToGp{this};
  GpToHpParameters gpToHp{this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "eastward_wind",
    "hydrostatic_pressure_levels",
    "northward_wind",
    "unbalanced_pressure_levels_minus_one"});}
};

}  // namespace saber
