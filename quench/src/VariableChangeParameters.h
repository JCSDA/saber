/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/VariableChangeParametersBase.h"

namespace quench {

// -------------------------------------------------------------------------------------------------
/// VariableChange parameters class

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, VariableChangeParametersBase)
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
