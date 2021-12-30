/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCETRAINING_H_
#define SABER_OOPS_ERRORCOVARIANCETRAINING_H_

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/BUMP.h"
#include "saber/oops/GSI.h"
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

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Background state parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// Ensemble parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};

  /// Ensemble perturbations parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};

  /// Ensemble base parameters
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};

  /// Background error covariance model
  oops::OptionalParameter<CovarianceParameters_> backgroundError{"background error", this};

  /// Input variables
  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  /// BUMP training parameters
  oops::OptionalParameter<BUMP_Parameters<MODEL>> bumpParams{"bump", this};

  /// GSI training parameters
  oops::OptionalParameter<GSI_Parameters> gsiParams{"gsi", this};
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
  typedef GSI<MODEL>                                         GSI_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceTraining";}
  explicit ErrorCovarianceTraining(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }

  virtual ~ErrorCovarianceTraining() {}

  int execute(const eckit::Configuration & fullConfig) const {
    util::Timer timer(classname(), "write");

    // Deserialize parameters
    ErrorCovarianceTrainingParameters_ params;
    params.validateAndDeserialize(fullConfig);

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
    const boost::optional<eckit::LocalConfiguration> &ensembleConfig = params.ensemble.value();
    if (ensembleConfig != boost::none) {
      // Ensemble of state, compute perturbation using the mean
      oops::Log::info() << "Ensemble of state, compute perturbation using the mean" << std::endl;
      ens1.reset(new Ensemble_(*ensembleConfig, xx, xx, resol, inputVars));
    } else {
      const boost::optional<eckit::LocalConfiguration>
        &ensemblePertConfig = params.ensemblePert.value();
      if (ensemblePertConfig != boost::none) {
        const boost::optional<eckit::LocalConfiguration>
          &ensembleBaseConfig = params.ensembleBase.value();
        if (ensembleBaseConfig != boost::none) {
          // Increment ensemble from difference of two state ensembles
          oops::Log::info() << "Increment ensemble from difference of two state ensembles"
                            << std::endl;
          ens1.reset(new Ensemble_(resol, inputVars, *ensembleBaseConfig, *ensemblePertConfig));
        } else {
          // Increment ensemble from increments on disk
          oops::Log::info() << "Increment ensemble from increments on disk" << std::endl;
          ens1.reset(new Ensemble_(resol, inputVars, *ensemblePertConfig));
        }
      }
    }
    if (ens1) {
      for (size_t ie = 0; ie < ens1->size(); ++ie) {
        (*ens1)[ie].toAtlas();
      }
    }

    // Setup ensemble 2
    EnsemblePtr_ ens2 = NULL;
    const boost::optional<CovarianceParameters_> &covarParams = params.backgroundError.value();
    if (covarParams != boost::none) {
      // Covariance matrix
      const CovarianceParametersBase_ &covarParamsBase = (*covarParams).covarianceParameters;
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                       covarParamsBase, resol, inputVars, xx, xx));

      // Randomize ensemble
      int ens2_ne = covarParamsBase.randomizationSize.value();
      ens2.reset(new Ensemble_(resol, inputVars, time, ens2_ne));
      for (int ie = 0; ie < ens2_ne; ++ie) {
        oops::Log::info() << "Generate randomized ensemble member " << ie+1 << " / "
                          << ens2_ne << std::endl;
        Increment_ incr(resol, inputVars, time);
        Bmat->randomize(incr);
        (*ens2)[ie] = incr;
        (*ens2)[ie].toAtlas();
      }
    }

    // Select SABER library training

    // BUMP
    const boost::optional<BUMP_Parameters<MODEL>> &bumpParams = params.bumpParams.value();
    if (bumpParams != boost::none) {
      // Do training
      BUMP_ bump(resol, inputVars, *bumpParams, ens1, ens2);

      // Write training parameters
      bump.write();

      // Apply training operators
      bump.apply();
    }

    // GSI
    const boost::optional<GSI_Parameters> &gsiParams = params.gsiParams.value();
    if (gsiParams != boost::none) {
      // Do training
      GSI_ gsi(resol, inputVars, *gsiParams);
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
