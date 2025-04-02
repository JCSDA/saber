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
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberBlockParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SaberBlockParametersBase, Parameters)

 public:
  // REQUIRED
  // SABER block name
  oops::RequiredParameter<std::string> saberBlockName{"saber block name", this};

  // PARAMETERS WITH DEFAULTS
  // Variables to update with left inverse for background and first guess
  oops::Parameter<oops::Variables> inverseVars{"state variables to inverse", oops::Variables(),
    this};

  // Flag to skip inverse application (for ensemble perturbations, background and first guess)
  oops::Parameter<bool> skipInverse{"skip inverse", false, this};

  // Flag to skip inverse test
  oops::Parameter<bool> skipInverseTest{"skip inverse test", false, this};

  // Flag to run the left inverse instead of the adjoint.
  oops::Parameter<bool> filterMode{"filter mode", false, this};

  // Fields metadata (e.g. geographical mask name, vertical coordinate field name)
  oops::Parameter<eckit::LocalConfiguration> fieldsMetaData{"fields metadata",
    eckit::LocalConfiguration(), this};

  // OPTIONAL
  // Active variables
  oops::OptionalParameter<oops::Variables> activeVars{"active variables", this};

  // Adjoint tolerance
  oops::OptionalParameter<double> adjointTolerance{"adjoint tolerance", this};

  // Calibration of block parameters
  oops::OptionalParameter<eckit::LocalConfiguration> calibrationParams{"calibration", this};

  // Ensemble transform parameters for the Ensemble block
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleTransform{"ensemble transform", this};

  // Tolerance for inner inverse test (U Uinv (U x) == (U x))
  oops::OptionalParameter<double> innerInverseTolerance{"inner inverse tolerance", this};

  // Inner variables to compare in outer inverse test, default is all inner active variables.
  oops::OptionalParameter<oops::Variables> innerVariables{"inner variables to compare", this};

  // Localization parameters for the Ensemble block
  oops::OptionalParameter<eckit::LocalConfiguration> localization{"localization", this};

  // Tolerance for outer inverse test (Uinv U (Uinv x) == (Uinv x))
  oops::OptionalParameter<double> outerInverseTolerance{"outer inverse tolerance", this};

  // Outer variables to compare in inner inverse test, default is all outer active variables.
  oops::OptionalParameter<oops::Variables> outerVariables{"outer variables to compare", this};

  // Read block parameters
  oops::OptionalParameter<eckit::LocalConfiguration> readParams{"read", this};

  // Tolerance for square-root test (U U^t x) == B x)
  oops::OptionalParameter<double> sqrtTolerance{"square-root tolerance", this};

  // Force calling write() even in "read" mode
  oops::Parameter<bool> forceWrite{"force write", false, this};

  // METHODS
  // Find out whether calibration is needed
  bool doCalibration() const;

  // Find out whether read is needed
  bool doRead() const;

  // VIRTUAL METHODS
  // Mandatory active variables
  virtual oops::Variables mandatoryActiveVars() const = 0;

  // Mandatory active inner variables, must be a subset of mandatoryActiveVars()
  // Used by utilities getUnionOfInnerActiveAndOuterVars() and getInnerOnlyVars()
  virtual oops::Variables activeInnerVars(const oops::Variables & outerVars) const {
    return oops::Variables();
  }

  // Mandatory active outer variables, must be a subset of mandatoryActiveVars()
  // Can be used to define outer variables to allocate when randomizing
  virtual oops::Variables activeOuterVars(const oops::Variables & outerVars) const {
    return oops::Variables();
  }

  virtual const oops::Variables mandatoryStateVars() const {
    return oops::Variables();
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber
