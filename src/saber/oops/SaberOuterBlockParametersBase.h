/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberOuterBlockParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SaberOuterBlockParametersBase, Parameters)
 public:
  oops::RequiredParameter<std::string> saberBlockName{"saber block name", this};
  oops::OptionalParameter<oops::Variables> inputVars{"input variables", this};
  oops::OptionalParameter<oops::Variables> activeVars{"active variables", this};
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputFields{"input fields", this};

  oops::Parameter<oops::Variables> outputVars{"output variables", oops::Variables(), this};
};

// -----------------------------------------------------------------------------

}  // namespace saber
