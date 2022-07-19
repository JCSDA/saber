/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SABERBLOCKTEST_H_
#define SABER_OOPS_SABERBLOCKTEST_H_

#include <memory>
#include <string>
#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTestParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockTestParameters, oops::ApplicationParameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_    StateParameters_;
  typedef SaberBlockParametersWrapper<MODEL>          SaberBlockParametersWrapper_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Tested variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// SABER blocks
  oops::RequiredParameter<std::vector<SaberBlockParametersWrapper<MODEL>>>
    saberBlocks{"saber blocks", this};

  /// Adjoint test tolerance
  oops::Parameter<double> adjointTolerance{"adjoint test tolerance", 1.0e-12, this};

  /// Inverse test tolerance
  oops::Parameter<double> inverseTolerance{"inverse test tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef SaberBlockBase<MODEL>                           SaberBlockBase_;
  typedef SaberBlockParametersWrapper<MODEL>              SaberBlockParametersWrapper_;
  typedef typename boost::ptr_vector<SaberBlockBase_>     SaberBlockVec_;
  typedef typename SaberBlockVec_::iterator               iter_;
  typedef typename SaberBlockVec_::const_iterator         icst_;
  typedef typename SaberBlockVec_::const_reverse_iterator ircst_;
  typedef oops::State<MODEL>                              State_;
  typedef SaberBlockTestParameters<MODEL>                 SaberBlockTestParameters_;

 public:
  static const std::string classname() {return "saber::SaberBlockTest";}
  explicit SaberBlockTest(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}

  virtual ~SaberBlockTest() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    SaberBlockTestParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables vars(params.variables);

    // Setup background state
    const State_ xx(geom, params.background);

    // Create SABER blocks
    std::unique_ptr<SaberBlockBase_> saberCentralBlock_;
    SaberBlockVec_ saberBlocks_;
    for (const SaberBlockParametersWrapper_ & saberBlockParamWrapper :
         params.saberBlocks.value()) {
      const SaberBlockParametersBase & saberBlockParams =
        saberBlockParamWrapper.saberBlockParameters;
      if (saberBlockParams.saberCentralBlock.value()) {
        if (saberCentralBlock_ || (saberBlocks_.size() != 0)) {
          ABORT("Central block should be the first block, only one allowed!");
        } else {
          saberCentralBlock_.reset(SaberBlockFactory<MODEL>::create(geom, saberBlockParams, xx,
            xx));
        }
      } else {
        saberBlocks_.push_back(SaberBlockFactory<MODEL>::create(geom, saberBlockParams, xx, xx));
      }
    }

    // Create test increments
    Increment_ dx1(geom, vars, xx.validTime());
    Increment_ dx2(geom, vars, xx.validTime());
    Increment_ dx1save(geom, vars, xx.validTime());
    Increment_ dx2save(geom, vars, xx.validTime());

    // Adjoint test for central block
    if (saberCentralBlock_) {
      // Generate random increments and save them
      dx1.random();
      dx2.random();
      dx1save = dx1;
      dx2save = dx2;

      // Apply central block
      saberCentralBlock_->multiply(dx1.fieldSet());
      saberCentralBlock_->multiply(dx2.fieldSet());

      // ATLAS fieldset to Increment_
      dx1.synchronizeFields();
      dx2.synchronizeFields();

      // Compute adjoint test
      const double dp1 = dx1.dot_product_with(dx2save);
      const double dp2 = dx2.dot_product_with(dx1save);
      oops::Log::test() << "Adjoint test for central block " << saberCentralBlock_->name()
        << std::endl;
      ASSERT(0.5*abs(dp1-dp2)/(dp1+dp2) < params.adjointTolerance.value());
    }

    // Adjoint test for other blocks
    for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
      // Generate random increments and save them
      dx1.random();
      dx2.random();
      dx1save = dx1;
      dx2save = dx2;

      // Apply non-central blocks
      it->multiply(dx1.fieldSet());
      it->multiplyAD(dx2.fieldSet());

      // ATLAS fieldset to Increment_
      dx1.synchronizeFields();
      dx2.synchronizeFields();

      // Compute adjoint test
      const double dp1 = dx1.dot_product_with(dx2save);
      const double dp2 = dx2.dot_product_with(dx1save);
      oops::Log::test() << "Adjoint test for block " << it->name() << std::endl;
      ASSERT(0.5*abs(dp1-dp2)/(dp1+dp2) < params.adjointTolerance.value());
    }

    // Inverse test for central block
    if (saberCentralBlock_) {
      if (saberCentralBlock_->iterativeInverse()) {
        oops::Log::test() << "Inverse test for central block " << saberCentralBlock_->name()
          << ": not implemented for iterative inverse" << std::endl;
      } else {
        // Generate random increments and save them
        dx1.random();
        dx1save = dx1;

        // Apply central block
        saberCentralBlock_->multiply(dx1.fieldSet());
        saberCentralBlock_->inverseMultiply(dx1.fieldSet());

        // ATLAS fieldset to Increment_
        dx1.synchronizeFields();

        // Compute inverse test
        dx1 -= dx1save;
        const double dp1 = dx1.norm()/dx1save.norm();
        oops::Log::test() << "Inverse test for central block " << saberCentralBlock_->name()
          << std::endl;
        ASSERT(dp1 < params.inverseTolerance.value());
      }
    }

    // Inverse test for other blocks
    for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
      if (it->iterativeInverse()) {
        oops::Log::test() << "Inverse test for non-central block " << it->name()
          << ": not implemented for iterative inverse" << std::endl;
      } else {
      // Generate random increments and save them
        dx1.random();
        dx2.random();
        dx1save = dx1;
        dx2save = dx2;

        // Apply central block
        it->multiply(dx1.fieldSet());
        it->inverseMultiply(dx1.fieldSet());
        it->multiplyAD(dx2.fieldSet());
        it->inverseMultiplyAD(dx2.fieldSet());

        // ATLAS fieldset to Increment_
        dx1.synchronizeFields();
        dx2.synchronizeFields();

        // Compute adjoint test
        dx1 -= dx1save;
        dx2 -= dx2save;
        const double dp1 = dx1.norm()/dx1save.norm();
        const double dp2 = dx2.norm()/dx2save.norm();
        oops::Log::test() << "Inverse test for block " << it->name() << std::endl;
        ASSERT(dp1 < params.inverseTolerance.value());
        ASSERT(dp2 < params.inverseTolerance.value());
      }
    }

    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::SaberBlockTest<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKTEST_H_
