/*
 * (C) Copyright 2025- UCAR
 * (C) Copyright 2025- NOAA/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "oops/coupled/TraitCoupled.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/ErrorCovarianceParameters.h"

namespace saber {

/// @brief  Parameters class for interpolation configurations for coupled error covariance
template <typename MODEL1, typename MODEL2>
class InterpolationForCoupledErrorCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationForCoupledErrorCovarianceParameters, Parameters)
 public:
  oops::RequiredParameter<eckit::LocalConfiguration>
    interp1{MODEL1::name().c_str(), this};
  oops::RequiredParameter<eckit::LocalConfiguration>
    interp2{MODEL2::name().c_str(), this};
};

// -------------------------------------------------------------------------------------------------
/// Parameters class for block-diagonal SABER coupled error covariance
template <typename MODEL1, typename MODEL2>
class CoupledErrorCovarianceParameters : public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(CoupledErrorCovarianceParameters, ErrorCovarianceParametersBase)

 public:
  oops::RequiredParameter<ErrorCovarianceParameters>
    errorCov1{MODEL1::name().c_str(), this};
  oops::RequiredParameter<ErrorCovarianceParameters>
    errorCov2{MODEL2::name().c_str(), this};

  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    commonOuterBlocks{"common outer blocks", this};
  oops::OptionalParameter<InterpolationForCoupledErrorCovarianceParameters<MODEL1, MODEL2>>
    interp{"interpolation", this};
};

// -----------------------------------------------------------------------------

}  // namespace saber
