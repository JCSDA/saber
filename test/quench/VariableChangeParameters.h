/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/Printable.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// \brief Parameters passed to the VariableChange class.

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, VariableChangeParametersBase)
 public:
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
