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

#include "saber/blocks/BlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class GaussUVToGPParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPParameters, BlockParametersBase)

 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const override;

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override;

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars{{outerVars["geostrophic_pressure_levels_minus_one"]}};
    return vars;
  }
};

// -----------------------------------------------------------------------------

class GpToHpReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GpToHpReadParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::RequiredParameter<int> covariance_nlat{"number of covariance latitude rings", this};
  oops::Parameter<int> gp_regression_bins{"gp regression bins", "gP regression bins", 18, this};
};

// -----------------------------------------------------------------------------

class GpToHpCalibrationReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GpToHpCalibrationReadParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::RequiredParameter<int> covariance_nlat{"number of covariance latitude rings", this};
  oops::Parameter<int> gp_regression_bins{"gp regression bins", "gP regression bins", 18, this};
};

// -----------------------------------------------------------------------------

class GpToHpParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GpToHpParameters, BlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<GpToHpReadParameters> readParams{"read", this};

  oops::OptionalParameter<GpToHpCalibrationReadParameters>
    calibrationReadParams{"calibration read", this};

  oops::Variables mandatoryActiveVars() const override;

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({"air_pressure_levels"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override;

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["hydrostatic_pressure_levels"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class GpToHpm1Parameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GpToHpm1Parameters, BlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<GpToHpReadParameters> readParams{"read", this};

  oops::OptionalParameter<GpToHpCalibrationReadParameters>
    calibrationReadParams{"calibration read", this};

  oops::Variables mandatoryActiveVars() const override;

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override;

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["hydrostatic_pressure_levels_minus_one"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class HydrostaticPressureParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticPressureParameters, BlockParametersBase)

 public:
  GaussUVToGPParameters gaussUVToGp{this};
  GpToHpParameters gpToHp{this};
  oops::Variables mandatoryActiveVars() const override;

  // same as in GpToHp parameters since this block is used in HydrostaticPressure
  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({"air_pressure_levels"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override;

  // activeOuterVars() is not needed in this super-block.
  // It would have contained "hydrostatic_pressure_levels".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const;
};

// -----------------------------------------------------------------------------

class HydrostaticPressureMinusOneParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticPressureMinusOneParameters, BlockParametersBase)

 public:
  GaussUVToGPParameters gaussUVToGp{this};
  GpToHpm1Parameters gpToHp{this};
  oops::Variables mandatoryActiveVars() const override;

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override;

  // activeOuterVars() is not needed in this super-block.
  // It would have contained "hydrostatic_pressure_levels_minus_one".

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const;
};

}  // namespace saber
