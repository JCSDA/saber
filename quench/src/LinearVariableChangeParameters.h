/*
 * (C) Copyright 2021 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LinearVariableChangeParametersBase.h"

#include "oops/base/Variables.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// LinearVariableChange parameters class

class LinearVariableChangeParameters : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, LinearVariableChangeParametersBase)
 public:
  // Variables (input/output)
  oops::OptionalParameter<oops::Variables> variables{"variables", this};

  // ATLAS file (multiplicative factor)
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFile{"atlas file", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
