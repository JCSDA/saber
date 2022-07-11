/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SABERBLOCKPARAMETERSBASE_H_
#define SABER_OOPS_SABERBLOCKPARAMETERSBASE_H_

#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberBlockParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SaberBlockParametersBase, Parameters)
 public:
  /// Common parameters.
  oops::RequiredParameter<std::string> saberBlockName{"saber block name", this};
  oops::Parameter<bool> saberCentralBlock{"saber central block", false, this};
  oops::Parameter<bool> iterativeInverse{"iterative inverse", false, this};
  oops::RequiredParameter<oops::Variables> inputVars{"input variables", this};
  oops::RequiredParameter<oops::Variables> outputVars{"output variables", this};
  oops::OptionalParameter<oops::Variables> activeVars{"active variables", this};
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKPARAMETERSBASE_H_
