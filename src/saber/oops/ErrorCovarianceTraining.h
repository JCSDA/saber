/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCETRAINING_H_
#define SABER_OOPS_ERRORCOVARIANCETRAINING_H_

#include <memory>
#include <string>

#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/instantiateCovarFactory.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

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

  /// Ensemble parameters
  oops::OptionalParameter<IncrementEnsembleFromStatesParameters_> ensemble{"ensemble", this};

  /// Ensemble perturbations parameters
  oops::OptionalParameter<IncrementEnsembleParameters_> ensemblePert{"ensemble pert", this};

  /// Ensemble base parameters
  oops::OptionalParameter<StateEnsembleParameters_> ensembleBase{"ensemble base", this};
  /// Ensemble state parameters for the ensemble pairs that would be subtracted from the base
  /// ensemble
  oops::OptionalParameter<StateEnsembleParameters_> ensemblePairs{"ensemble pairs", this};

  /// Background error covariance model
  oops::OptionalParameter<CovarianceParameters_> backgroundError{"background error", this};

  /// Input variables
  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  /// BUMP training parameters
  oops::OptionalParameter<BUMP_Parameters<MODEL>> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ErrorCovarianceTraining : public oops::Application {
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef oops::ModelSpaceCovarianceBase<MODEL>              CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                     CovarianceFactory_;
  typedef oops::ModelSpaceCovarianceParametersBase<MODEL>    CovarianceParametersBase_;
  typedef oops::Geometry<MODEL>                              Geometry_;
  typedef oops::Increment<MODEL>                             Increment_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::IncrementEnsemble<MODEL>                     Ensemble_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>>    EnsemblePtr_;
  typedef ErrorCovarianceTrainingParameters<MODEL>           ErrorCovarianceTrainingParameters_;
  typedef BUMP<MODEL>                                        BUMP_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceTraining";}
  explicit ErrorCovarianceTraining(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }

  virtual ~ErrorCovarianceTraining() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    ErrorCovarianceTrainingParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    //  Setup resolution
    const Geometry_ resol(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables inputVars(params.inputVariables);

    // Setup background state
    const State_ xx(resol, params.background);

    // Setup time
    const util::DateTime time = xx.validTime();

    // Setup ensemble 1
    EnsemblePtr_ ens1 = NULL;
    if (params.ensemble.value() != boost::none) {
      // Ensemble of state, compute perturbation using the mean
      oops::Log::info() << "Ensemble of state, compute perturbation using the mean" << std::endl;
      ens1.reset(new Ensemble_(*params.ensemble.value(), xx, xx, resol, inputVars));
    } else if (params.ensemblePert.value() != boost::none) {
      // Increment ensemble from increments on disk
      oops::Log::info() << "Increment ensemble from increments on disk" << std::endl;
      ens1.reset(new Ensemble_(resol, inputVars, *params.ensemblePert.value()));
    } else if ((params.ensembleBase.value() != boost::none) &&
               (params.ensemblePairs.value() != boost::none)) {
      // Increment ensemble from difference of two state ensembles
       oops::Log::info() << "Increment ensemble from difference of two state ensembles"
                         << std::endl;
       ens1.reset(new Ensemble_(resol, inputVars, *params.ensembleBase.value(),
                                                  *params.ensemblePairs.value()));
    }

    // Setup ensemble 2
    EnsemblePtr_ ens2 = NULL;
    const boost::optional<CovarianceParameters_> &covarParams = params.backgroundError.value();
    if (covarParams != boost::none) {
      // Covariance matrix
      const CovarianceParametersBase_ &covarParamsBase = (*covarParams).covarianceParameters;
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                       resol, inputVars, covarParamsBase, xx, xx));

      // Randomize ensemble
      int ens2_ne = covarParamsBase.randomizationSize.value();
      ens2.reset(new Ensemble_(resol, inputVars, time, ens2_ne));
      for (int ie = 0; ie < ens2_ne; ++ie) {
        oops::Log::info() << "Generate randomized ensemble member " << ie+1 << " / "
                          << ens2_ne << std::endl;
        Increment_ incr(resol, inputVars, time);
        Bmat->randomize(incr);
        (*ens2)[ie] = incr;
      }
    }

    // Select SABER library training

    // BUMP
    const boost::optional<BUMP_Parameters<MODEL>> &bumpParams = params.bumpParams.value();
    if (bumpParams != boost::none) {
      // Do training
      BUMP_ bump(resol, inputVars, *bumpParams, xx, xx, ens1, ens2);

      // Write training parameters
      bump.write();

      // Apply training operators
      bump.apply();
    }

    return 0;
  }

 private:
  std::string appname() const {
    return "saber::ErrorCovarianceTraining<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCETRAINING_H_
