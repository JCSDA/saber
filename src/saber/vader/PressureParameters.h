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

class GaussUVToGPParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "eastward_wind",
    "geostrophic_pressure_levels_minus_one",
    "northward_wind"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"eastward_wind",
                          "northward_wind"});
    const int modelLevels = outerVars.getLevels("geostrophic_pressure_levels_minus_one");
    vars.addMetaData("eastward_wind", "levels", modelLevels);
    vars.addMetaData("northward_wind", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"geostrophic_pressure_levels_minus_one"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

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
  oops::RequiredParameter<GpToHpCovarianceParameters>
    gptohpcovarianceparams{"covariance data", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "geostrophic_pressure_levels_minus_one",
    "hydrostatic_pressure_levels",
    "unbalanced_pressure_levels_minus_one"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"geostrophic_pressure_levels_minus_one",
                          "unbalanced_pressure_levels_minus_one"});
    const int modelLevels = outerVars.getLevels("hydrostatic_pressure_levels") - 1;
    vars.addMetaData("geostrophic_pressure_levels_minus_one", "levels", modelLevels);
    vars.addMetaData("unbalanced_pressure_levels_minus_one", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"hydrostatic_pressure_levels"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

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

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"eastward_wind",
                          "northward_wind",
                          "unbalanced_pressure_levels_minus_one"});
    const int modelLevels = outerVars.getLevels("hydrostatic_pressure_levels") - 1;
    vars.addMetaData("eastward_wind", "levels", modelLevels);
    vars.addMetaData("northward_wind", "levels", modelLevels);
    vars.addMetaData("unbalanced_pressure_levels_minus_one", "levels", modelLevels);
    return vars;
  }

  // activeOuterVars() is not needed in this super-block.
  // It would have contained "hydrostatic_pressure_levels".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
    oops::Variables tempVars({"geostrophic_pressure_levels_minus_one"});
    if (outerVars.has("geostrophic_pressure_levels_minus_one")) {
      throw eckit::UserError("geostrophic_pressure_levels_minus_one is a "
                             "temporary variable of mo_hydrostatic_pressure "
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels = outerVars.getLevels("hydrostatic_pressure_levels") - 1;
    tempVars.addMetaData("geostrophic_pressure_levels_minus_one", "levels", modelLevels);
    return tempVars;
  }
};

}  // namespace saber
