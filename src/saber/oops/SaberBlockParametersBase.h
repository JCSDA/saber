/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberBlockParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SaberBlockParametersBase, Parameters)

 public:
  // SABER block name
  oops::RequiredParameter<std::string> saberBlockName{"saber block name", this};

  // Active variables [optional]
  oops::OptionalParameter<oops::Variables> activeVars{"active variables", this};

  // Read block parameters [optional]
  oops::OptionalParameter<eckit::LocalConfiguration> readParams{"read", this};

  // Calibration of block parameters [optional]
  oops::OptionalParameter<eckit::LocalConfiguration> calibrationParams{"calibration", this};

  // Localization parameter for the Ensemble block
  oops::OptionalParameter<eckit::LocalConfiguration> localization{"localization", this};

  // Adjoint tolerance
  oops::OptionalParameter<double> adjointTolerance{"adjoint tolerance", this};

  // Flag to skip inverse test
  oops::Parameter<bool> skipInverseTest{"skip inverse test", false, this};

  // Tolerance for inner inverse test (U Uinv (U x) == (U x))
  oops::OptionalParameter<double> innerInverseTolerance{"inner inverse tolerance", this};

  // Tolerance for outer inverse test (Uinv U (Uinv x) == (Uinv x))
  oops::OptionalParameter<double> outerInverseTolerance{"outer inverse tolerance", this};

  // Inner variables to compare in outer inverse test, default is all inner active variables.
  oops::OptionalParameter<oops::Variables> innerVariables{"inner variables to compare", this};

  // Outer variables to compare in inner inverse test, default is all outer active variables.
  oops::OptionalParameter<oops::Variables> outerVariables{"outer variables to compare", this};

  // Find out whether calibration is needed
  bool doCalibration() const;

  // Find out whether read is needed
  bool doRead() const;

  // Mandatory active variables
  virtual oops::Variables mandatoryActiveVars() const = 0;
};

// -----------------------------------------------------------------------------

}  // namespace saber
