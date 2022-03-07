/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SABERBLOCKTEST_H_
#define SABER_OOPS_SABERBLOCKTEST_H_

#include <memory>
#include <string>
#include <vector>

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

#include "saber/oops/instantiateSaberBlockFactory.h"
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
  typedef typename oops::Geometry<MODEL>::Parameters_  GeometryParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Input variables
  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  /// Dirac location/variables parameters
  oops::RequiredParameter<eckit::LocalConfiguration> dirac{"dirac", this};

  /// SABER blocks
  oops::RequiredParameter<std::vector<SaberBlockParametersWrapper<MODEL>>>
      saberBlocks{"saber blocks", this};

  /// Output parameters
  oops::RequiredParameter<eckit::LocalConfiguration> output{"output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::State<MODEL>                              State_;
  typedef SaberBlockBase<MODEL>                           SaberBlockBase_;
  typedef SaberBlockParametersWrapper<MODEL>              SaberBlockParametersWrapper_;
  typedef SaberBlockTestParameters<MODEL> SaberBlockTestParameters_;
  typedef typename boost::ptr_vector<SaberBlockBase_>     SaberBlockVec_;
  typedef typename SaberBlockVec_::iterator               iter_;
  typedef typename SaberBlockVec_::const_iterator         icst_;
  typedef typename SaberBlockVec_::const_reverse_iterator ircst_;

 public:
  static const std::string classname() {return "saber::SaberBlockTest";}
  explicit SaberBlockTest(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateSaberBlockFactory<MODEL>();
  }

  virtual ~SaberBlockTest() {}

  int execute(const eckit::Configuration & fullConfig) const {
    // Deserialize parameters
    SaberBlockTestParameters_ params;
    params.validateAndDeserialize(fullConfig);
    oops::Log::info() << "fullConfig: " << fullConfig << std::endl;

    //  Setup resolution
    const Geometry_ resol(params.geometry.value(), this->getComm());

    // Setup variables
    const oops::Variables inputVars(params.inputVariables.value());
    oops::Log::info() << "Variables: " << inputVars << std::endl;

    // Dummy date
    const util::DateTime date(1977, 5, 25, 0, 0, 0);
    oops::Log::info() << "Date: " << date << std::endl;

    // Setup state
    State_ xx(resol, inputVars, date);

    // Setup increments
    Increment_ dx(resol, inputVars, date);

    // Initial dirac increment
    dx.dirac(params.dirac.value());
    oops::Log::test() << "Initial increment: " << dx << std::endl;

    // Check input/output variables consistency
    oops::Variables vars_in(inputVars);
    oops::Variables vars_out;
    for (const SaberBlockParametersWrapper_ & saberBlockParamWrapper :
         boost::adaptors::reverse(params.saberBlocks.value())) {
      const SaberBlockParametersBase & saberBlockParams =
        saberBlockParamWrapper.saberBlockParameters;
      vars_out = saberBlockParams.outputVars.value();
      if (!(vars_in == vars_out)) {
        const boost::optional<std::string> &saberBlockName =
          saberBlockParams.saberBlockName.value();
        if (saberBlockName != boost::none) {
          oops::Log::error() << "In SABER block " << *saberBlockName << std::endl;
        } else {
          oops::Log::error() << "In SABER block with no name" << std::endl;
        }
        oops::Log::error() << "  Input variables:  " << vars_in << std::endl;
        oops::Log::error() << "  Output variables: " << vars_out << std::endl;
        ABORT("  Sequence of blocks is not consistent (wrong variables)");
      }
      vars_in = saberBlockParams.inputVars.value();
    }

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
          saberCentralBlock_.reset(SaberBlockFactory<MODEL>::create(resol, saberBlockParams, xx,
            xx));
        }
      } else {
        saberBlocks_.push_back(SaberBlockFactory<MODEL>::create(resol, saberBlockParams, xx, xx));
      }
    }

    // Increment_ to ATLAS fieldset
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());

    // K_1^T K_2^T .. K_N^T
    for (ircst_ it = saberBlocks_.rbegin(); it != saberBlocks_.rend(); ++it) {
      it->multiplyAD(atlasFieldSet.get());
    }

    // Central block multiplication
    if (saberCentralBlock_) {
      saberCentralBlock_->multiply(atlasFieldSet.get());
    }

    // K_N K_N-1 ... K_1
    for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
      it->multiply(atlasFieldSet.get());
    }

    // Write increment
    oops::Log::test() << "Final increment: " << dx << std::endl;
    dx.write(params.output.value());
    return 0;
  }

 private:
  std::string appname() const {
    return "saber::SaberBlockTest<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKTEST_H_
