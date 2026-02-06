/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {

// -------------------------------------------------------------------------------------------------

class ErrorCovarianceParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParametersBase, oops::Parameters)

 public:
  // Adjoint test
  oops::Parameter<bool> adjointTest{"adjoint test", false, this};
  oops::Parameter<double> adjointTolerance{"adjoint tolerance", 1.0e-12, this};

  // Inverse test
  oops::Parameter<bool> inverseTest{"inverse test", false, this};
  oops::Parameter<double> inverseTolerance{"inverse tolerance", 1.0e-12, this};

  // Square-root test
  oops::Parameter<bool> sqrtTest{"square-root test", false, this};
  oops::Parameter<double> sqrtTolerance{"square-root tolerance", 1.0e-12, this};
};

// -------------------------------------------------------------------------------------------------

class ErrorCovarianceParameters : public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParameters,
                           ErrorCovarianceParametersBase)

 public:
  /// ModelSpaceCovarianceBase class parameters

  // Covariance model
  oops::OptionalParameter<std::string> covarianceModel{"covariance model", this};

  // Randomization size
  oops::OptionalParameter<size_t> randomizationSize{"randomization size", this};

  // Inverse parameters
  oops::Parameter<bool> fullInverse{"full inverse", false, this};
  oops::Parameter<int> fullInverseIterations{"full inverse iterations", 10, this};
  oops::Parameter<double> fullInverseAccuracy{"full inverse accuracy", 1.0e-3, this};

  // Extra linear variable change
  oops::OptionalParameter<eckit::LocalConfiguration> variableChange{"linear variable change", this};

  /// SABER ErrorCovariance-specific parameters

  // Option to change resolution of the background to the increment geometry
  oops::Parameter<bool> changeBackgroundResolution{"change background resolution",
                        false, this};

  /// Block chain parameters
  oops::ConfigurationParameter blockChainParams{this};
};

// -----------------------------------------------------------------------------

}  // namespace saber
