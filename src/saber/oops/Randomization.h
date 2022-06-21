/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_RANDOMIZATION_H_
#define SABER_OOPS_RANDOMIZATION_H_

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
    Increment_ dx(geom, vars, xx.validTime());
    for (size_t jm = 0; jm < covarParams.randomizationSize.value(); ++jm) {
      // Generate pertubation
      Bmat->randomize(dx);
      oops::Log::test() << "Member " << jm << ": " << dx << std::endl;

      // Add mean state
      State_ xp(xx);
      xp += dx;

      // Write perturbation
      StateWriterParameters_ outParams = params.output;
      outParams.write.setMember(jm+1);
      xp.write(outParams.write);
    }

    return 0;
  }

 private:
  std::string appname() const {
    return "saber::Randomization<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_RANDOMIZATION_H_
