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

class GaussUVToGPParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "eastward_wind",
    "geostrophic_pressure_levels_minus_one",
    "northward_wind"}});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["geostrophic_pressure_levels_minus_one"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back(oops::Variable{"eastward_wind", conf});
    vars.push_back(oops::Variable{"northward_wind", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars{{outerVars["geostrophic_pressure_levels_minus_one"]}};
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
    std::vector<std::string>{
    "geostrophic_pressure_levels_minus_one",
    "hydrostatic_pressure_levels",
    "unbalanced_pressure_levels_minus_one"}});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["hydrostatic_pressure_levels"].getLevels() - 1;
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back(oops::Variable{"geostrophic_pressure_levels_minus_one", conf});
    vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["hydrostatic_pressure_levels"]});
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
    std::vector<std::string>{"eastward_wind",
    "hydrostatic_pressure_levels",
    "northward_wind",
    "unbalanced_pressure_levels_minus_one"}});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["hydrostatic_pressure_levels"].getLevels() - 1;
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back(oops::Variable{"eastward_wind", conf});
    vars.push_back(oops::Variable{"northward_wind", conf});
    vars.push_back(oops::Variable{"unbalanced_pressure_levels_minus_one", conf});
    return vars;
  }

  // activeOuterVars() is not needed in this super-block.
  // It would have contained "hydrostatic_pressure_levels".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
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
};

}  // namespace saber
