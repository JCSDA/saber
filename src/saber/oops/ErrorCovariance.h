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

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

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
  oops::RequiredParameter<std::vector<SaberBlockParametersWrapper<MODEL>>>
      saberBlocks{"saber blocks", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef SaberBlockBase<MODEL>                           SaberBlockBase_;
  typedef SaberBlockParametersWrapper<MODEL>              SaberBlockParametersWrapper_;
  typedef typename boost::ptr_vector<SaberBlockBase_>     SaberBlockVec_;
  typedef typename SaberBlockVec_::iterator               iter_;
  typedef typename SaberBlockVec_::const_iterator         icst_;
  typedef typename SaberBlockVec_::const_reverse_iterator ircst_;
  typedef oops::State<MODEL>                              State_;

 public:
  typedef ErrorCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const Parameters_ &,
                  const State_ &, const State_ &);
  virtual ~ErrorCovariance();

  // Required by iterative inverses for central blocks
  void multiply(const Increment_ &, Increment_ &) const;

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<SaberBlockBase_> saberCentralBlock_;
  SaberBlockVec_ saberBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & resol,
                                        const oops::Variables & inputVars,
                                        const Parameters_ & params,
                                        const State_ & xb, const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(resol, params, xb, fg), saberCentralBlock_(),
    saberBlocks_()
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Check input/output variables consistency
  oops::Variables vars_in(inputVars);
  oops::Variables vars_out;

  for (const SaberBlockParametersWrapper_ & saberBlockParamWrapper :
       boost::adaptors::reverse(params.saberBlocks.value())) {
    const SaberBlockParametersBase & saberBlockParams = saberBlockParamWrapper.saberBlockParameters;
    vars_out = saberBlockParams.outputVars.value();
    if (!(vars_in == vars_out)) {
      const boost::optional<std::string> &saberBlockName = saberBlockParams.saberBlockName.value();
      if (saberBlockName != boost::none) {
        oops::Log::error() << "In SABER block \'" << *saberBlockName << "\'" << std::endl;
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
  for (const SaberBlockParametersWrapper_ & saberBlockParamWrapper :
       params.saberBlocks.value()) {
    const SaberBlockParametersBase & saberBlockParams = saberBlockParamWrapper.saberBlockParameters;
    if (saberBlockParams.saberCentralBlock.value()) {
      if (saberCentralBlock_ || (saberBlocks_.size() != 0)) {
        ABORT("Central block should be the first block, only one allowed!");
      } else {
        saberCentralBlock_.reset(SaberBlockFactory<MODEL>::create(resol, saberBlockParams, xb,
          fg));
      }
    } else {
      saberBlocks_.push_back(SaberBlockFactory<MODEL>::create(resol, saberBlockParams, xb, fg));
    }
  }

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

  // Randomization done flag
  bool randDone(false);

  // Central block randomization
  if (saberCentralBlock_) {
    saberCentralBlock_->randomize(dx.fieldSet());
    randDone = true;
  }

  // K_N K_N-1 ... K_1
  for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
    if (!randDone) {
      it->randomize(dx.fieldSet());
      randDone = true;
    } else {
      it->multiply(dx.fieldSet());
    }
  }

  // ATLAS fieldset to Increment_
  dx.synchronizeFieldsAD();

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

  // K_1^T K_2^T .. K_N^T
  for (ircst_ it = saberBlocks_.rbegin(); it != saberBlocks_.rend(); ++it) {
    it->multiplyAD(dxo.fieldSet());
  }

  // Central block multiplication
  if (saberCentralBlock_) {
    saberCentralBlock_->multiply(dxo.fieldSet());
  }

  // K_N K_N-1 ... K_1
  for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
    it->multiply(dxo.fieldSet());
  }

  // ATLAS fieldset to Increment_
  dxo.synchronizeFieldsAD();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                               Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Copy input
  dxo = dxi;

  // K_1^{-1} K_2^{-1} .. K_N^{-1}
  for (ircst_ it = saberBlocks_.rbegin(); it != saberBlocks_.rend(); ++it) {
    it->inverseMultiply(dxo.fieldSet());
  }

  // Synchronization flag
  bool syncNeeded(true);

  // Central block inverse multiplication
  if (saberCentralBlock_) {
    if (saberCentralBlock_->iterativeInverse()) {
      if (saberBlocks_.size() > 0) {
        // ATLAS fieldset to Increment_
        dxo.synchronizeFieldsAD();
      } else {
        syncNeeded = false;
      }

      // Temporary increment
      Increment_ dxtmp(dxo);

      // Iterative inverse
      oops::IdentityMatrix<Increment_> Id;
      dxo.zero();
      GMRESR(dxo, dxtmp, *this, Id, 10, 1.0e-3);
    } else {
      // Block-specific inverse
      saberCentralBlock_->inverseMultiply(dxo.fieldSet());
    }
  }

  // K_N^T^{-1} K_N-1^T^{-1} ... K_1^T^{-1}
  for (icst_ it = saberBlocks_.begin(); it != saberBlocks_.end(); ++it) {
    it->inverseMultiplyAD(dxo.fieldSet());
  }

  // ATLAS fieldset to Increment_
  if (syncNeeded) {
    dxo.synchronizeFieldsAD();
  }

  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::multiply(const Increment_ & dxi,
                                      Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Copy input
  dxo = dxi;

  // Central block multiplication
  if (saberCentralBlock_) {
    saberCentralBlock_->multiply(dxo.fieldSet());
  } else {
    ABORT("iterative inverse for central blocks only");
  }

  // ATLAS fieldset to Increment_
  dxo.synchronizeFieldsAD();

  oops::Log::trace() << "ErrorCovariance<MODEL>::multiply done" << std::endl;
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
