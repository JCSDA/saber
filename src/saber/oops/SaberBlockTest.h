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

#include "boost/range/adaptors.hpp"
#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Random.h"

#include "saber/oops/ReadInputFields.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

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

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Tested variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// SABER blocks
  oops::OptionalParameter<SaberCentralBlockParametersWrapper>
    saberCentralBlock{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocks{"saber outer blocks", this};

  /// Adjoint test tolerance
  oops::Parameter<double> adjointTolerance{"adjoint test tolerance", 1.0e-12, this};

  /// Inverse test tolerance
  oops::Parameter<double> inverseTolerance{"inverse test tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef typename boost::ptr_vector<SaberOuterBlockBase>      SaberOuterBlockVec_;
  typedef typename SaberOuterBlockVec_::iterator               iter_;
  typedef typename SaberOuterBlockVec_::const_iterator         icst_;
  typedef typename SaberOuterBlockVec_::const_reverse_iterator ircst_;
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

    // Saber blocks
    std::unique_ptr<SaberCentralBlockBase> saberCentralBlock_;
    SaberOuterBlockVec_ saberOuterBlocks_;

    // Local copy of background and first guess
    State_ xbLocal(xx);
    State_ fgLocal(xx);

    // Initial and output variables and geometry
    oops::Variables outputVars(incVars);
    atlas::FunctionSpace outputFunctionSpace = geom.functionSpace();
    atlas::FieldSet outputExtraFields = geom.extraFields();

    // Build outer blocks successively
    const boost::optional<std::vector<SaberOuterBlockParametersWrapper>> &saberOuterBlocks = params.saberOuterBlocks.value();
    if (saberOuterBlocks != boost::none) {
      // Loop in reverse order
      for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
        boost::adaptors::reverse(*saberOuterBlocks)) {
        // Get outer block parameters
        const SaberOuterBlockParametersBase & saberOuterBlockParams = saberOuterBlockParamWrapper.saberOuterBlockParameters;
  
        // Local configuration to add parameters
        eckit::LocalConfiguration outerConf;
        saberOuterBlockParams.serialize(outerConf);
        outerConf.set("output variables", outputVars.variables());
  
        // Define active variables
        oops::Variables activeVars;
        const boost::optional<oops::Variables> &optionalActiveVars = saberOuterBlockParams.activeVars.value();
        if (optionalActiveVars != boost::none) {
           // Active variables specified
           activeVars = *optionalActiveVars;
        } else {
           // No active variables specified, assuming they are the same as output variables
           activeVars = outputVars;
           outerConf.set("active variables", activeVars.variables());
        }

        // Read input fields (on model increment geometry)
        std::vector<atlas::FieldSet> fsetVec = readInputFields(
          geom,
          activeVars,
          xx.validTime(),
          saberOuterBlockParams.inputFields.value());
  
        // Create outer block
        oops::Log::info() << "Creating outer block: " << saberOuterBlockParams.saberBlockName.value() << std::endl;
        saberOuterBlocks_.push_back(SaberOuterBlockFactory::create(geom.getComm(),
                                    outputFunctionSpace,
                                    outputExtraFields,
                                    geom.variableSizes(activeVars),
                                    outerConf,
                                    xbLocal.fieldSet(),
                                    fgLocal.fieldSet(),
                                    fsetVec));
  
        // Apply calibration inverse on xb and fg
        saberOuterBlocks_.back().calibrationInverseMultiply(xbLocal.fieldSet());
        saberOuterBlocks_.back().calibrationInverseMultiply(fgLocal.fieldSet());
  
        // Access input geometry and variables of the current block
        const oops::Variables inputVars = saberOuterBlocks_.back().inputVars();
        const atlas::FunctionSpace inputFunctionSpace = saberOuterBlocks_.back().inputFunctionSpace();
        const atlas::FieldSet inputExtraFields = saberOuterBlocks_.back().inputExtraFields();

        // Check that active variables are present in either input or output variables, or both
        for (const auto & var : activeVars.variables()) {
          ASSERT(inputVars.has(var) || outputVars.has(var));
        }

        // Adjoint test

        // Variables sizes
        std::vector<size_t> variableSizes = geom.variableSizes(inputVars);

        // Create random input FieldSet
        atlas::FieldSet inputFset = createRandomFieldSet(inputFunctionSpace,
                                                         variableSizes,
                                                         inputVars);
/*
        // Save input and output Fieldsets

        // Apply adjoint block
        saberBlock_->multiplyAD();

      // 
      Increment_ dx2TLAD(dx1TLAD);
      dx2TLAD.random();
      Increment_ dx2TLADsave(dx2TLAD, true);
  
      // Apply forward block
      saberBlock_->multiply(dx2TLAD.fieldSet());

      // ATLAS fieldset to Increment_
      dx1TLAD.synchronizeFieldsAD();
      dx2TLAD.synchronizeFields();

      // Compute adjoint test
      const double dp1 = dx1TLAD.dot_product_with(dx2TLADsave);
      const double dp2 = dx2TLAD.dot_product_with(dx1TLADsave);

      oops::Log::test() << "Adjoint test for block " << saberCentralBlock_->name() <<
                           ": y^t (Ax) = " << dp1 <<
                           ": x^t (Ay) = " << dp2 << std::endl;
      ASSERT(0.5*abs(dp1-dp2)/(dp1+dp2) < params.adjointTolerance.value());

  */
        // Update output geometry and variables for the next block
        outputFunctionSpace = inputFunctionSpace;
        outputExtraFields = inputExtraFields;
        outputVars = inputVars;
      }
    }

    const boost::optional<SaberCentralBlockParametersWrapper> &saberCentralBlock = params.saberCentralBlock.value();
    if (saberCentralBlock != boost::none) {
      // Get central block parameters
      const SaberCentralBlockParametersBase & saberCentralBlockParams = saberCentralBlock->saberCentralBlockParameters;
    
      // Define input/output variables
      oops::Variables inoutVars = outputVars;
    
      // Local configuration to add parameters
      eckit::LocalConfiguration centralConf;
      saberCentralBlockParams.serialize(centralConf);
      centralConf.set("inout variables", inoutVars.variables());
    
      // Define active variables
      oops::Variables activeVars;
      const boost::optional<oops::Variables> &optionalActiveVars = saberCentralBlockParams.activeVars.value();
      if (optionalActiveVars != boost::none) {
         // Active variables specified
         activeVars = *optionalActiveVars;
      } else {
         // No active variables specified, assuming they are the same as output variables
         activeVars = inoutVars;
         centralConf.set("active variables", activeVars.variables());
      }
    
      // Read input fields (on model increment geometry)
      std::vector<atlas::FieldSet> fsetVec = readInputFields(
        geom,
        activeVars,
        xx.validTime(),
        saberCentralBlockParams.inputFields.value());
    
      // Create central block
      saberCentralBlock_.reset(SaberCentralBlockFactory::create(geom.getComm(),
                               outputFunctionSpace,
                               outputExtraFields,
                               geom.variableSizes(activeVars),
                               centralConf,
                               xbLocal.fieldSet(),
                               fgLocal.fieldSet(),
                               fsetVec));
    
      // Check that active variables are present in input/output variables
      for (const auto & var : activeVars.variables()) {
        ASSERT(inoutVars.has(var));
      }
    }

    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::SaberBlockTest<" + MODEL::name() + ">";
  }

  // -----------------------------------------------------------------------------
  // TODO: should be in OOPS?
  atlas::FieldSet createRandomFieldSet(const atlas::FunctionSpace & functionSpace,
                                       const std::vector<size_t> & variableSizes,
                                       const oops::Variables & vars) {
  
    // Get ghost points
    atlas::Field ghost = functionSpace.ghost();
    auto ghostView = atlas::array::make_view<int, 1>(ghost);
  
    // Create FieldSet
    atlas::FieldSet fset;
    for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
      // Create field
      atlas::Field field = functionSpace.createField<double>(
        atlas::option::name(vars.variables()[jvar])
        | atlas::option::levels(variableSizes[jvar]));
  
      // Get field owned size
      size_t n = 0;
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (ghostView(jnode) == 0) ++n;
          }
        }
      }

      // Generate random vector
      util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);
  
      // Populate with random numbers
      n = 0;
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (ghostView(jnode) == 0) {
              view(jnode, jlevel) = rand_vec[n];
              ++n;
            } else {
              view(jnode, jlevel) = 0.0;
            }
          }
        }
      }       
  
      // Add field
      fset.add(field);
    }
  
    // Halo exchange
    fset.haloExchange();
  
    // Return FieldSet
    return fset;
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKTEST_H_
