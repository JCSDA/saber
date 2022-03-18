/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_COVARIANCEPARAMS_H_
#define QUENCH_COVARIANCEPARAMS_H_

#include "oops/base/ModelSpaceCovarianceParametersBase.h"

#include "quench/LinearVarChange.h"
#include "quench/TraitsFwd.h"

namespace quench {
// -----------------------------------------------------------------------------
/// \brief Parameters passed to the CovarianceParams class.

class CovarianceParams : public oops::ModelSpaceCovarianceParametersBase<quench::Traits>{
  OOPS_CONCRETE_PARAMETERS(CovarianceParams,
                           oops::ModelSpaceCovarianceParametersBase<quench::Traits>)

 public:
};
// -----------------------------------------------------------------------------

}  // namespace quench
#endif  // QUENCH_COVARIANCEPARAMS_H_
