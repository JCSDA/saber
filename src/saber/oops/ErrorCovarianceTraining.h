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

  /// Randomized ensemble output
  oops::OptionalParameter<eckit::LocalConfiguration>
    randomizedEnsembleOutput{"randomized ensemble output", this};

  /// Geometry parameters for ensemble 2
  oops::OptionalParameter<GeometryParameters_> geometry2{"lowres geometry", this};

  /// Ensemble 2 parameters
  oops::OptionalParameter<IncrementEnsembleFromStatesParameters_> ensemble2{"lowres ensemble",
    this};

  /// Ensemble 2 perturbations parameters
  oops::OptionalParameter<IncrementEnsembleParameters_> ensemble2Pert{"lowres ensemble pert",
    this};

  /// Ensemble 2 base parameters
  oops::OptionalParameter<StateEnsembleParameters_> ensemble2Base{"lowres ensemble base",
    this};

  /// Ensemble 2 state parameters for the ensemble pairs that would be subtracted from
  /// the base ensemble
  oops::OptionalParameter<StateEnsembleParameters_> ensemble2Pairs{"lowres ensemble pairs",
    this};

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

    // Setup geometry
    const Geometry_ geom1(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables inputVars(params.inputVariables);

    // Setup background state
    const State_ xx(geom1, params.background);

    // Setup time
    const util::DateTime time = xx.validTime();

    // Setup ensemble 1
    EnsemblePtr_ ens1 = NULL;
    const boost::optional<IncrementEnsembleFromStatesParameters_>
      &ensemble = params.ensemble.value();
    const boost::optional<IncrementEnsembleParameters_> &ensemblePert = params.ensemblePert.value();
    const boost::optional<StateEnsembleParameters_> &ensembleBase = params.ensembleBase.value();
    const boost::optional<StateEnsembleParameters_> &ensemblePairs = params.ensemblePairs.value();
    if (ensemble != boost::none) {
      // Ensemble of states, perturbation using the mean
      oops::Log::info() << "Ensemble of states, perturbation using the mean" << std::endl;
      ens1.reset(new Ensemble_(*ensemble, xx, xx, geom1, inputVars));
    } else if (ensemblePert) {
      // Increment ensemble from increments on disk
      oops::Log::info() << "Increment ensemble from increments on disk" << std::endl;
      ens1.reset(new Ensemble_(geom1, inputVars, *ensemblePert));
    } else if ((ensembleBase != boost::none) &&
               (ensemblePairs != boost::none)) {
      // Increment ensemble from difference of two states
       oops::Log::info() << "Increment ensemble from difference of two states"
                         << std::endl;
       ens1.reset(new Ensemble_(geom1, inputVars, *ensembleBase, *ensemblePairs));
    }

    // Setup ensemble 2 geometry pointer
    const Geometry_ * geom2 = &geom1;
    const boost::optional<GeometryParameters_> &geom2Params = params.geometry2.value();
    if (geom2Params != boost::none) {
      geom2 = new Geometry_(*geom2Params, geom1.getComm());
    }

    // Setup ensemble 2
    EnsemblePtr_ ens2 = NULL;
    const boost::optional<CovarianceParameters_>
      &backgroundError = params.backgroundError.value();
    const boost::optional<IncrementEnsembleFromStatesParameters_>
      &ensemble2 = params.ensemble2.value();
    const boost::optional<IncrementEnsembleParameters_>
      &ensemble2Pert = params.ensemble2Pert.value();
    const boost::optional<StateEnsembleParameters_> &ensemble2Base = params.ensemble2Base.value();
    const boost::optional<StateEnsembleParameters_> &ensemble2Pairs = params.ensemble2Pairs.value();

    if (backgroundError != boost::none) {
      // Covariance matrix
      const CovarianceParametersBase_ &covarParamsBase = (*backgroundError).covarianceParameters;
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                       *geom2, inputVars, covarParamsBase, xx, xx));

      // Randomize ensemble and remove mean
      const int ens2_ne = covarParamsBase.randomizationSize.value();
      ens2.reset(new Ensemble_(*geom2, inputVars, time, ens2_ne));
      Increment_ mean(*geom2, inputVars, time);
      mean.zero();
      for (int ie = 0; ie < ens2_ne; ++ie) {
        oops::Log::info() << "Generate randomized ensemble member " << ie+1 << " / "
                          << ens2_ne << std::endl;
        Increment_ incr(*geom2, inputVars, time);
        Bmat->randomize(incr);
        (*ens2)[ie] = incr;
        mean += incr;
      }
      const double rr = 1.0/static_cast<double>(ens2_ne);
      mean *= rr;
      for (int ie = 0; ie < ens2_ne; ++ie) {
        (*ens2)[ie] -= mean;
      }

      // Optionally write randomized ensemble
      const boost::optional<eckit::LocalConfiguration>
        &randomizedEnsembleOutput = params.randomizedEnsembleOutput.value();
      if (randomizedEnsembleOutput != boost::none) {
        ens2->write(*randomizedEnsembleOutput);
      }
    } else if ((ensemble2 != boost::none) || (ensemble2Pert != boost::none)
     || ((ensemble2Base != boost::none) && (ensemble2Pairs != boost::none))) {
      // Setup low resolution background state
      const State_ xx2(*geom2, xx);

      // Low resolution ensemble
      if (ensemble2 != boost::none) {
        // Low resolution ensemble of states, perturbation using the mean
        oops::Log::info() << "Low resolution ensemble of states, perturbation using the mean"
          << std::endl;
        ens2.reset(new Ensemble_(*ensemble2, xx2, xx2, *geom2, inputVars));
      } else if (ensemble2Pert != boost::none) {
        // Low resolution increment ensemble from increments on disk
        oops::Log::info() << "Low resolution increment ensemble from increments on disk"
          << std::endl;
        ens2.reset(new Ensemble_(*geom2, inputVars, *ensemble2Pert));
      } else if ((ensemble2Base != boost::none) && (ensemble2Pairs != boost::none)) {
        // Low resolution increment ensemble from difference of two states
         oops::Log::info() << "Low resolution increment ensemble from difference of two states"
           << std::endl;
         ens2.reset(new Ensemble_(*geom2, inputVars, *ensemble2Base, *ensemble2Pairs));
      }
    }

    // Select SABER library training

    // BUMP
    const boost::optional<BUMP_Parameters<MODEL>> &bumpParams = params.bumpParams.value();
    if (bumpParams != boost::none) {
      // Do training
      BUMP_ bump(geom1, *geom2, inputVars, *bumpParams, xx, xx, ens1, ens2);
    }

    // Delete pointer
    if (geom2Params != boost::none) {
      delete geom2;
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
