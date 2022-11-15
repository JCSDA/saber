/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/ModelSpaceCovarianceParametersBase.h"

#include "quench/LinearVariableChange.h"
#include "quench/TraitsFwd.h"

namespace quench {
// -----------------------------------------------------------------------------
/// \brief Parameters passed to the Covariance class.

class CovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<quench::Traits>{
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<quench::Traits>)

 public:
};
// -----------------------------------------------------------------------------

}  // namespace quench
