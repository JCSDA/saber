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

#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bump/BUMPParameters.h"

#include "saber/bump/lib/BUMP.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class InputParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InputParameters, oops::Parameters)
  typedef typename oops::Increment<MODEL>::ReadParameters_ ReadParameters_;

 public:
  /// Parameter name.
  oops::RequiredParameter<std::string> param{"parameter", this};

  /// Component index
  oops::Parameter<int> component{"component", 1, this};

  /// Parameters used for reading Increment.
  oops::RequiredParameter<ReadParameters_> file{"file", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class OutputParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OutputParameters, oops::Parameters)
  typedef typename oops::Increment<MODEL>::WriteParameters_ WriteParameters_;

 public:
  /// Parameter name.
  oops::RequiredParameter<std::string> param{"parameter", this};

  /// Component index
  oops::Parameter<int> component{"component", 1, this};

  /// Parameters used for writing Increment.
  oops::RequiredParameter<WriteParameters_> file{"file", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ErrorCovarianceTrainingParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceTrainingParameters, oops::ApplicationParameters)

 public:
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef typename oops::Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_           StateParameters_;
  typedef oops::StateEnsembleParameters<MODEL>               StateEnsembleParameters_;
  typedef oops::IncrementEnsembleFromStatesParameters<MODEL> IncrementEnsembleFromStatesParameters_;
  typedef oops::IncrementEnsembleParameters<MODEL>           IncrementEnsembleParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Background state parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// Background error covariance model
  oops::RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ErrorCovarianceTraining : public oops::Application {
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef oops::ModelSpaceCovarianceBase<MODEL>              CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                     CovarianceFactory_;
  typedef oops::ModelSpaceCovarianceParametersBase<MODEL>    CovarianceParametersBase_;
  typedef oops::StateEnsembleParameters<MODEL>               StateEnsembleParameters_;
  typedef oops::IncrementEnsembleFromStatesParameters<MODEL> IncrementEnsembleFromStatesParameters_;
  typedef oops::IncrementEnsembleParameters<MODEL>           IncrementEnsembleParameters_;
  typedef typename oops::Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef oops::Geometry<MODEL>                              Geometry_;
  typedef oops::Increment<MODEL>                             Increment_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::IncrementEnsemble<MODEL>                     Ensemble_;
  typedef std::shared_ptr<Ensemble_>                         EnsemblePtr_;
  typedef ErrorCovarianceTrainingParameters<MODEL>           ErrorCovarianceTrainingParameters_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceTraining";}
  explicit ErrorCovarianceTraining(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }

  virtual ~ErrorCovarianceTraining() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Get number of MPI tasks and OpenMP threads
    std::string mpi(std::to_string(this->getComm().size()));
    std::string omp("1");
#ifdef _OPENMP
    # pragma omp parallel
    {
      omp = std::to_string(omp_get_num_threads());
    }
#endif

    // Deserialize parameters
    ErrorCovarianceTrainingParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom1(params.geometry, this->getComm());

    // Setup background state
    const State_ xx(geom1, params.background);

    // Create covariance
    const CovarianceParametersBase_ &covarParams =
        params.backgroundError.value().covarianceParameters;
    eckit::LocalConfiguration covarConf(covarParams.toConfiguration());
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
        xx.geometry(), xx.variables(), covarConf, xx, xx));

    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::ErrorCovarianceTraining<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber
