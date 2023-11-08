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

#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {

// -------------------------------------------------------------------------------------------------

class DualResCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DualResCalibrationParameters, oops::Parameters)
 public:
  // Geometry
  oops::OptionalParameter<eckit::LocalConfiguration> geometry{"geometry", this};

  // Ensemble
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<MODEL>)

 public:
  // Central and outer blocks
  oops::RequiredParameter<SaberCentralBlockParametersWrapper>
    saberCentralBlockParams{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocksParams{"saber outer blocks", this};

  // Time covariance mode (by default duplicated multivariate)
  // Options: univariate, duplicated multivariate.
  oops::Parameter<std::string> timeCovariance{"time covariance", "multivariate duplicated",
                                              this};

  // Ensemble
  oops::Parameter<bool> iterativeEnsembleLoading{"iterative ensemble loading", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};

  // Dual resolution calibration
  oops::OptionalParameter<DualResCalibrationParameters> dualResParams{
    "dual resolution calibration", this};

  // Output ensemble
  oops::OptionalParameter<eckit::LocalConfiguration> outputEnsemble{"output ensemble", this};

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

// -----------------------------------------------------------------------------

}  // namespace saber
