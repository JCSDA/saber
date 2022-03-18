/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_LINEARVARCHANGEPARAMS_H_
#define QUENCH_LINEARVARCHANGEPARAMS_H_

#include <ostream>
#include <string>

#include "oops/base/LinearVariableChangeParametersBase.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// No change of linear variable parameters

class LinearVarChangeParams : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVarChangeParams, LinearVariableChangeParametersBase)
 public:
  // Linear variable change parameters would go here.
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
#endif  // QUENCH_LINEARVARCHANGEPARAMS_H_
