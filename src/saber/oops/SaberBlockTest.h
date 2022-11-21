/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/range/adaptors.hpp>

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/ReadInputFields.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockWrapper.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTestParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockTestParameters, oops::ApplicationParameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_    StateParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Tested variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// SABER blocks
  oops::OptionalParameter<SaberCentralBlockWrapperParameters<MODEL>>
    saberCentralBlock{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocks{"saber outer blocks", this};

  /// Adjoint test tolerance
  oops::Parameter<double> adjointTolerance{"adjoint test tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef typename boost::ptr_vector<SaberOuterBlockBase>      SaberOuterBlockVec_;
  typedef oops::State<MODEL>                                   State_;

 public:
  static const std::string classname() {return "saber::SaberBlockTest";}
  explicit SaberBlockTest(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}

  virtual ~SaberBlockTest() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    SaberBlockTestParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables incVars(params.variables);

    // Setup background state
    const State_ xx(geom, params.background);

    // SABER blocks
    std::unique_ptr<SaberCentralBlockWrapper<MODEL>> saberCentralBlock_;
    SaberOuterBlockVec_ saberOuterBlocks_;

    // Local copy of background and first guess
    State_ xbLocal(xx);
    State_ fgLocal(xx);

    // Initial outer geometry and variables
    std::vector<std::reference_wrapper<const oops::GeometryData>> outerGeometryData_;
    outerGeometryData_.push_back(geom.generic());
    oops::Variables outerVars(incVars);

    // Build outer blocks successively
    const boost::optional<std::vector<SaberOuterBlockParametersWrapper>> &saberOuterBlocks =
      params.saberOuterBlocks.value();
    if (saberOuterBlocks != boost::none) {
      // Loop in reverse order
      for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
        boost::adaptors::reverse(*saberOuterBlocks)) {
        // Get outer block parameters
        const SaberBlockParametersBase & saberOuterBlockParams =
          saberOuterBlockParamWrapper.saberOuterBlockParameters;

        // Get active variables
        oops::Variables activeVars =
          saberOuterBlockParams.activeVars.value().get_value_or(outerVars);

        // Read input fields (on model increment geometry)
        std::vector<eckit::LocalConfiguration> inputFields;
        inputFields = saberOuterBlockParams.inputFields.value().get_value_or(inputFields);
        std::vector<atlas::FieldSet> fsetVec = readInputFields(geom,
                                                               activeVars,
                                                               xx.validTime(),
                                                               inputFields);

        // Create outer block
        oops::Log::info() << "Info     : Creating outer block: "
                          << saberOuterBlockParams.saberBlockName.value() << std::endl;
        saberOuterBlocks_.push_back(SaberOuterBlockFactory::create(
                                    outerGeometryData_.back().get(),
                                    geom.variableSizes(activeVars),
                                    outerVars,
                                    saberOuterBlockParams,
                                    xbLocal.fieldSet(),
                                    fgLocal.fieldSet(),
                                    fsetVec));

        // Access inner geometry and variables
        const oops::GeometryData & innerGeometryData = saberOuterBlocks_.back().innerGeometryData();
        const oops::Variables innerVars = saberOuterBlocks_.back().innerVars();

        // Check that active variables are present in either inner or outer variables, or both
        for (const auto & var : activeVars.variables()) {
          ASSERT(innerVars.has(var) || outerVars.has(var));
        }

        // Adjoint test

        // Variables sizes
        std::vector<size_t> innerVariableSizes = geom.variableSizes(innerVars);

        // Create random inner FieldSet
        atlas::FieldSet innerFset = util::createRandomFieldSet(innerGeometryData.functionSpace(),
                                                         innerVariableSizes,
                                                         innerVars);

        // Copy inner FieldSet
        atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

        // Variables sizes
        std::vector<size_t> outerVariableSizes = geom.variableSizes(outerVars);

        // Create random outer FieldSet
        atlas::FieldSet outerFset =
          util::createRandomFieldSet(outerGeometryData_.back().get().functionSpace(),
                               outerVariableSizes,
                               outerVars);

        // Copy outer FieldSet
        atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

        // Apply forward and adjoint multiplication
        saberOuterBlocks_.back().multiply(innerFset);
        saberOuterBlocks_.back().multiplyAD(outerFset);

        // Compute adjoint test
        const double dp1 = util::dotProductFieldSets(innerFset, outerFsetSave, activeVars,
                                                     geom.getComm());
        const double dp2 = util::dotProductFieldSets(outerFset, innerFsetSave, activeVars,
                                                     geom.getComm());
        oops::Log::info() << "Info     : Adjoint test for outer block "
                          << saberOuterBlockParams.saberBlockName.value()
                          << ": y^t (Ax) = " << dp1 << ": x^t (A^t y) = " << dp2 << std::endl;
        ASSERT(abs(dp1) > 0.0);
        ASSERT(abs(dp2) > 0.0);
        oops::Log::test() << "Adjoint test for outer block "
                          << saberOuterBlockParams.saberBlockName.value();
        if (0.5*abs(dp1-dp2)/(dp1+dp2) < params.adjointTolerance.value()) {
          oops::Log::test() << " passed" << std::endl;
        } else {
          oops::Log::test() << " failed" << std::endl;
          ABORT("Adjoint test failure");
        }

        // Update outer geometry and variables for the next block
        outerGeometryData_.push_back(innerGeometryData);
        outerVars = innerVars;
      }
    }

    // Get central block wrapper parameters
    const boost::optional<SaberCentralBlockWrapperParameters<MODEL>> &saberCentralBlockParams =
      params.saberCentralBlock.value();

    if (saberCentralBlockParams != boost::none) {
      // Create central block wrapper
      saberCentralBlock_.reset(new SaberCentralBlockWrapper<MODEL>(geom,
                               outerGeometryData_.back().get(),
                               outerVars,
                               *saberCentralBlockParams,
                               xbLocal,
                               fgLocal,
                               params.adjointTolerance.value()));
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
