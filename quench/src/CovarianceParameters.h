/*
 * (C) Copyright 2021 UCAR
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/ModelSpaceCovarianceParametersBase.h"

#include "src/TraitsFwd.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Covariance parameters class

class CovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<Traits>{
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<Traits>)
};

// -----------------------------------------------------------------------------

}  // namespace quench
