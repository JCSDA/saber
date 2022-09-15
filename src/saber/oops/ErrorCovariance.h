/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCE_H_
#define SABER_OOPS_ERRORCOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/oops/ReadInputFields.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<MODEL>)
 public:
  oops::RequiredParameter<SaberCentralBlockParametersWrapper>
    saberCentralBlocks{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocks{"saber outer blocks", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef typename boost::ptr_vector<SaberOuterBlockBase>      SaberOuterBlockVec_;
  typedef typename SaberOuterBlockVec_::iterator               iter_;
  typedef typename SaberOuterBlockVec_::const_iterator         icst_;
  typedef typename SaberOuterBlockVec_::const_reverse_iterator ircst_;
  typedef oops::State<MODEL>                                   State_;

 public:
  typedef ErrorCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const Parameters_ &,
                  const State_ &, const State_ &);
  virtual ~ErrorCovariance();

  // Required by iterative inverse
  void multiply(const Increment_ & dxi, Increment_ & dxo) const {this->doMultiply(dxi, dxo);}

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<SaberCentralBlockBase> saberCentralBlock_;
  SaberOuterBlockVec_ saberOuterBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol,
                                        const oops::Variables & incVars,
                                        const Parameters_ & params,
                                        const State_ & xb,
                                        const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(resol, params, xb, fg), saberCentralBlock_(),
    saberOuterBlocks_()
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Local copy of background and first guess
  State_ xbLocal(xb);
  State_ fgLocal(fg);

  // Initial and output variables and geometry
  oops::Variables outputVars(incVars);
  atlas::FunctionSpace outputFunctionSpace = resol.functionSpace();
  atlas::FieldSet outputExtraFields = resol.extraFields();

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

      // Define input variables
      oops::Variables inputVars;
      const boost::optional<oops::Variables> &optionalInputVars = saberOuterBlockParams.inputVars.value();
      if (optionalInputVars != boost::none) {
         // Input variables specified
         inputVars = *optionalInputVars;
      } else {
         // No input variables specified, assuming they are the same as output variables
         inputVars = outputVars;
         outerConf.set("input variables", inputVars.variables());
      }

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

      // Check that active variables are present
      for (const auto & var : activeVars.variables()) {
        ASSERT(inputVars.has(var) || outputVars.has(var));
      }

      // Define input geometry (TODO: access it from the block, only copy right now)
      atlas::FunctionSpace inputFunctionSpace = outputFunctionSpace;
      atlas::FieldSet inputExtraFields = outputExtraFields;

      // Read input fields (on model increment geometry)
      std::vector<atlas::FieldSet> fsetVec = readInputFields(
        resol,
        inputVars,
        xb.validTime(),
        saberOuterBlockParams.inputFields.value());

      // Create outer block
      oops::Log::info() << "Creating outer block: " << saberOuterBlockParams.saberBlockName.value() << std::endl;
      saberOuterBlocks_.push_back(SaberOuterBlockFactory::create(resol.getComm(),
                                  inputFunctionSpace,
                                  inputExtraFields,
                                  resol.variableSizes(inputVars),
                                  outputFunctionSpace,
                                  outputExtraFields,
                                  resol.variableSizes(outputVars),
                                  outerConf,
                                  xbLocal.fieldSet(),
                                  fgLocal.fieldSet(),
                                  fsetVec));

      // Apply calibration inverse on xb and fg
      saberOuterBlocks_.back().calibrationInverseMultiply(xbLocal.fieldSet());
      saberOuterBlocks_.back().calibrationInverseMultiply(fgLocal.fieldSet());

      // Input variables and geometry of this block will be output variables and geometry of the next one
      outputVars = inputVars;
      outputFunctionSpace = inputFunctionSpace;
      outputExtraFields = inputExtraFields;
    }
  }

  // Get central block parameters
  const SaberCentralBlockParametersWrapper & saberCentralBlockParamWrapper = params.saberCentralBlocks.value();
  const SaberCentralBlockParametersBase & saberCentralBlockParams = saberCentralBlockParamWrapper.saberCentralBlockParameters;

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

  // Check that active variables are present
  for (const auto & var : activeVars.variables()) {
    ASSERT(inoutVars.has(var));
  }

  // Read input fields (on model increment geometry)
  std::vector<atlas::FieldSet> fsetVec = readInputFields(
    resol,
    inoutVars,
    xb.validTime(),
    saberCentralBlockParams.inputFields.value());

  // Create central block
  saberCentralBlock_.reset(SaberCentralBlockFactory::create(resol.getComm(),
                           outputFunctionSpace,
                           outputExtraFields,
                           resol.variableSizes(inoutVars),
                           centralConf,
                           xbLocal.fieldSet(),
                           fgLocal.fieldSet(),
                           fsetVec));
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doRandomize(Increment_ & dx) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  // Random output vector (necessary for some SABER blocks)
  dx.random();

  // Central block randomization
  saberCentralBlock_->randomize(dx.fieldSet());

  // Outer blocks forward multiplication
  for (ircst_ it = saberOuterBlocks_.rbegin(); it != saberOuterBlocks_.rend(); ++it) {
    it->multiply(dx.fieldSet());
  }

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment_ & dxi,
                                        Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  // Copy input
  dxo = dxi;

  // Outer blocks adjoint multiplication
  for (icst_ it = saberOuterBlocks_.begin(); it != saberOuterBlocks_.end(); ++it) {
    it->multiplyAD(dxo.fieldSet());
  }

  // Central block multiplication
  if (saberCentralBlock_) {
    saberCentralBlock_->multiply(dxo.fieldSet());
  }

  // Outer blocks forward multiplication
  for (ircst_ it = saberOuterBlocks_.rbegin(); it != saberOuterBlocks_.rend(); ++it) {
    it->multiply(dxo.fieldSet());
  }

  // ATLAS fieldset to Increment_
  dxo.synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                               Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCE_H_
