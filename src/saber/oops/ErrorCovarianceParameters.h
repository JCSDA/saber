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

class DualResolutionCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DualResolutionCalibrationParameters, oops::Parameters)
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
class OutputEnsembleParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OutputEnsembleParameters, oops::Parameters)
 public:
  // Output ensemble
  typename oops::Increment<MODEL>::WriteParameters_ writeParams{this};

  // Write first member only
  oops::Parameter<bool> firstMemberOnly{"first member only", false, this};
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

  // Ensemble
  oops::Parameter<bool> iterativeEnsembleLoading{"iterative ensemble loading", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};

  // Dual resolution calibration
  oops::OptionalParameter<DualResolutionCalibrationParameters> dualResolutionParams{
    "dual resolution calibration", this};

  // Output ensemble
  oops::OptionalParameter<OutputEnsembleParameters<MODEL>> outputEnsemble{"output ensemble", this};

  // Adjoint test
  oops::Parameter<bool> adjointTest{"adjoint test", false, this};
  oops::Parameter<double> adjointTolerance{"adjoint tolerance", 1.0e-12, this};

  // Inverse test
  oops::Parameter<bool> inverseTest{"inverse test", false, this};
  oops::Parameter<double> inverseTolerance{"inverse tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

}  // namespace saber
