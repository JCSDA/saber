/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/StateWriter.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/instantiateCovarFactory.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class RandomizationParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(RandomizationParameters, oops::ApplicationParameters)

 public:
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef typename oops::Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_           StateParameters_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::StateWriterParameters<State_>                StateWriterParameters_;
  typedef typename oops::Increment<MODEL>::WriteParameters_  WriteParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Randomized variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background state parameters.
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// Background error covariance model
  oops::RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};

  /// Where to write the output
  oops::RequiredParameter<StateWriterParameters_> output{"output", this};

  /// Where to write the output
  oops::OptionalParameter<WriteParameters_> outputIncrement{"output increment", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class Randomization : public oops::Application {
  typedef oops::ModelSpaceCovarianceBase<MODEL>              CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                     CovarianceFactory_;
  typedef oops::ModelSpaceCovarianceParametersBase<MODEL>    CovarianceParametersBase_;
  typedef oops::Geometry<MODEL>                              Geometry_;
  typedef oops::Increment<MODEL>                             Increment_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::StateWriterParameters<State_>                StateWriterParameters_;
  typedef RandomizationParameters<MODEL>                     RandomizationParameters_;
  typedef typename oops::Increment<MODEL>::WriteParameters_  WriteParameters_;

 public:
  static const std::string classname() {return "saber::Randomization";}
  explicit Randomization(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }

  virtual ~Randomization() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    RandomizationParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables vars(params.variables);

    // Setup background state
    const State_ xx(geom, params.background);

    // Setup B matrix
    const CovarianceParametersBase_ &covarParams =
        params.backgroundError.value().covarianceParameters;
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                          geom, vars, covarParams, xx, xx));

    // Generate and write perturbations
    std::vector<Increment_> ens;
    Increment_ mean(geom, vars, xx.validTime());
    mean.zero();
    const boost::optional<WriteParameters_> &inputIncrOutParams = params.outputIncrement.value();
    for (size_t jm = 0; jm < covarParams.randomizationSize.value(); ++jm) {
      // Generate pertubation
      Increment_ dx(geom, vars, xx.validTime());
      Bmat->randomize(dx);
      oops::Log::test() << "Member " << jm << ": " << dx << std::endl;

      // Update mean
      mean += dx;

      // Save perturbation
      ens.push_back(dx);
    }

    // Normalize mean
    const double norm = 1.0/static_cast<double>(covarParams.randomizationSize.value());
    mean *= norm;

    for (size_t jm = 0; jm < covarParams.randomizationSize.value(); ++jm) {
      // Remove mean
      ens[jm] -= mean;

      // Write perturbation
      if (inputIncrOutParams != boost::none) {
        WriteParameters_ incrOutParams = *inputIncrOutParams;
        incrOutParams.setMember(jm+1);
        ens[jm].write(incrOutParams);
      }

      // Add mean state
      State_ xp(xx);
      xp += ens[jm];

      // Write perturbed state
      StateWriterParameters_ outParams = params.output;
      outParams.write.setMember(jm+1);
      xp.write(outParams.write);
    }
    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::Randomization<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber
