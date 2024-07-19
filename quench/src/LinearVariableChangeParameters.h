/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LinearVariableChangeParametersBase.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// \brief Parameters passed to the LinearVariableChange class.

class LinearVariableChangeParameters : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, LinearVariableChangeParametersBase)
 public:
  /// Variables (input/output)
  oops::OptionalParameter<oops::Variables> variables{"variables", this};

  /// ATLAS file
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFile{"atlas file", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
