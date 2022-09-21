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
#include <vector>

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
#include "saber/oops/ReadInputFields.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

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
  WriteParameters_ incwrite{this};
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
  oops::RequiredParameter<oops::Variables> inputVars{"input variables", this};

  /// Input fields 1
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputFields{"input fields", this};

  /// Input fields 2
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>>
    inputFields2{"lowres input fields", this};

  /// BUMP training parameters
  oops::OptionalParameter<BUMP_Parameters> bumpParams{"bump", this};

  /// Output parameters
  oops::OptionalParameter<std::vector<OutputParameters<MODEL>>> output{"output", this};
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

    // Deserialize parameters
    ErrorCovarianceTrainingParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom1(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables inputVars(params.inputVars);

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

    // Get input fields for geometry 1
    std::vector<atlas::FieldSet> fsetVec1 = readInputFields(
      geom1,
      params.inputVars.value(),
      xx.validTime(),
      params.inputFields.value());

    // Get input fields for geometry 2
    std::vector<atlas::FieldSet> fsetVec2 = readInputFields(
      geom1,
      params.inputVars.value(),
      xx.validTime(),
      params.inputFields2.value());

    // Select SABER library training
    std::unique_ptr<BUMP> bump;

    // Ensemble sizes
    size_t ens1_ne = 0;
    size_t ens2_ne = 0;
    if (ens1) ens1_ne = ens1->size();
    if (ens2) ens2_ne = ens2->size();

    // BUMP
    const boost::optional<BUMP_Parameters> &bumpParams = params.bumpParams.value();
    if (bumpParams != boost::none) {
      // Constructor
      bump.reset(new BUMP(geom1.getComm(),
                          geom1.functionSpace(),
                          geom1.extraFields(),
                          geom1.variableSizes(inputVars),
                          inputVars,
                          *bumpParams,
                          fsetVec1,
                          geom2->functionSpace(),
                          geom2->extraFields(),
                          fsetVec2,
                          ens1_ne,
                          ens2_ne));

      // Add members of ensemble 1
      if (ens1) {
        oops::Log::info() << "--- Add members of ensemble 1" << std::endl;
        for (size_t ie = 0; ie < ens1_ne; ++ie) {
          oops::Log::info() << "      Member " << ie+1 << " / " << ens1_ne << std::endl;
          bump->addMember((*ens1)[ie].fieldSet(), ie, 1);
        }
      }

      // Add members of ensemble 2
      if (ens2) {
        oops::Log::info() << "--- Add members of ensemble 2" << std::endl;
        for (size_t ie = 0; ie < ens2_ne; ++ie) {
          oops::Log::info() << "      Member " << ie+1 << " / " << ens2_ne << std::endl;
          bump->addMember((*ens2)[ie].fieldSet(), ie, 2);
        }
      }

      // Check what needs to be updated
      const boost::optional<bool> &update_vbal_cov = bumpParams->update_vbal_cov.value();
      const boost::optional<bool> &update_var = bumpParams->update_var.value();
      const boost::optional<bool> &update_mom = bumpParams->update_mom.value();

      // Load ensemble members sequentially
      if (bump->memberConfig1().size() > 0) {
        ens1_ne = bump->memberConfig1().size();
        Increment_ dx1(geom1, inputVars, xx.validTime());

        for (size_t ie = 0; ie < ens1_ne; ++ie) {
          // Read member
          oops::Log::info() <<
        "-------------------------------------------------------------------" << std::endl;
          oops::Log::info() << "--- Load member " << ie+1 << " / " << ens1_ne << std::endl;
          dx1.read(bump->memberConfig1()[ie]);

          if (update_vbal_cov != boost::none) {
            if (*update_vbal_cov) {
              // Update vertical covariance
              bump->updateVbalCov(dx1.fieldSet(), ie);
            }
          }
          if (update_var != boost::none) {
            if (*update_var) {
              // Update variance
              bump->updateVar(dx1.fieldSet(), ie);
            }
          }
          if (update_mom != boost::none) {
            if (*update_mom) {
              // Update moments
              bump->updateMom(dx1.fieldSet(), ie, 1);
            }
          }
        }
      }
      if (bump->memberConfig2().size() > 0) {
        ens2_ne = bump->memberConfig2().size();
        Increment_ dx2(*geom2, inputVars, xx.validTime());

        for (size_t ie = 0; ie < ens2_ne; ++ie) {
          // Read member
          oops::Log::info() <<
          "-------------------------------------------------------------------" << std::endl;
          oops::Log::info() << "--- Load member " << ie+1 << " / " << ens2_ne << std::endl;
          dx2.read(bump->memberConfig2()[ie]);
          if (update_mom != boost::none) {
            if (*update_mom) {
              // Update moments
              bump->updateMom(dx2.fieldSet(), ie, 2);
            }
          }
        }
      }

      // Run drivers
      bump->runDrivers();

      // Partial deallocation
      bump->partialDealloc();

      // Apply operators
      const boost::optional<std::vector<eckit::LocalConfiguration>>
        &appConfs = bumpParams->appConfs.value();
      if (appConfs != boost::none) {
        oops::Log::info() <<
        "-------------------------------------------------------------------" << std::endl;
        oops::Log::info() << "--- Apply operators" << std::endl;
        if (appConfs->size() > 0) {
          for (const auto & appConf : *appConfs) {
            // Read input file
            eckit::LocalConfiguration inputConf(appConf, "input");
            oops::Log::info() << "       - Input file: " << inputConf << std::endl;
            Increment_ dx1(geom1, inputVars, xx.validTime());
            dx1.read(inputConf);

            // Apply BUMP operator
            std::vector<std::string> bumpOperators;
            appConf.get("bump operators", bumpOperators);
            for (const auto & bumpOperator : bumpOperators) {
              oops::Log::info() << "         Apply operator " << bumpOperator << std::endl;
              if (bumpOperator == "multiplyVbal") {
                bump->multiplyVbal(dx1.fieldSet());
              } else if (bumpOperator == "inverseMultiplyVbal") {
                bump->inverseMultiplyVbal(dx1.fieldSet());
              } else if (bumpOperator == "multiplyVbalAd") {
                bump->multiplyVbalAd(dx1.fieldSet());
              } else if (bumpOperator == "inverseMultiplyAd") {
                bump->inverseMultiplyVbalAd(dx1.fieldSet());
              } else if (bumpOperator == "multiplyStdDev") {
                bump->multiplyStdDev(dx1.fieldSet());
              } else if (bumpOperator == "inverseMultiplyStdDev") {
                bump->inverseMultiplyStdDev(dx1.fieldSet());
              } else if (bumpOperator == "multiplyNicas") {
                bump->multiplyNicas(dx1.fieldSet());
              } else {
                  ABORT("Wrong bump operator: " + bumpOperator);
              }
            }

            // ATLAS fieldset to Increment_
            dx1.synchronizeFields();

            // Write file
            eckit::LocalConfiguration outputConf(appConf, "output");
            oops::Log::info() << "         Output file: " << outputConf << std::endl;
            dx1.write(outputConf);
          }
        }
      }
    }

    // Write output parameters to file
    const boost::optional<std::vector<OutputParameters<MODEL>>> &output = params.output.value();
    if (output != boost::none) {
      for (const auto & outputParam : *output) {
        // Get parameter
        const std::string & param = outputParam.param;

        // Get component
        const int & component = outputParam.component;

        // BUMP output
        if (bumpParams != boost::none) {
          // Select geometry
          if (param == "loc_a_lr"
           || param == "loc_rh_lr"
           || param == "loc_rh1_lr"
           || param == "loc_rh2_lr"
           || param == "loc_rhc_lr"
           ||param == "loc_rv_lr"
           || param == "dirac_diag_loc_lr"
           || param == "nicas_norm_lr"
           || param == "dirac_nicas_lr"
           || param == "dirac_nicas_bens_lr") {
            // Get parameter
            Increment_ dx2(*geom2, inputVars, xx.validTime());
            dx2.zero(xx.validTime());
            bump->getParameter(param, component, 2, dx2.fieldSet());
            dx2.synchronizeFields();

            // Write parameter
            dx2.write(outputParam.incwrite);
            oops::Log::test() << "Norm of BUMP output parameter " << param << " - " << component
                              << ": " << dx2.norm() << std::endl;
          } else {
            // Get parameter
            Increment_ dx1(geom1, inputVars, xx.validTime());
            dx1.zero(xx.validTime());
            bump->getParameter(param, component, 1, dx1.fieldSet());
            dx1.synchronizeFields();

            // Write parameter
            dx1.write(outputParam.incwrite);
            oops::Log::test() << "Norm of BUMP output parameter " << param << " - " << component
                              << ": " << dx1.norm() << std::endl;
          }
        }
      }
    }

    // Delete pointer
    if (geom2Params != boost::none) {
      delete geom2;
    }

    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::ErrorCovarianceTraining<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCETRAINING_H_
