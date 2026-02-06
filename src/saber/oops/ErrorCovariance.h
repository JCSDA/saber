/*
 * (C) Copyright 2021-2023 UCAR
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef oops::Increment4D<MODEL>                             Increment4D_;
  typedef oops::State4D<MODEL>                                 State4D_;

 public:
  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const eckit::Configuration &,
                  const State4D_ &, const State4D_ &);
  virtual ~ErrorCovariance();

  void multiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {this->doMultiply(dxi, dxo);}

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;

  /// Chain of blocks (hybrid or ensemble or parametric)
  std::unique_ptr<SaberBlockChainBase> blockChain_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const eckit::Configuration & config,
                                        const State4D_ & xb,
                                        const State4D_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, config, xb, fg)
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "ErrorCovariance");

  // Deserialize parameters
  ErrorCovarianceParameters params;
  params.deserialize(config);

  // Local copy of background and first guess that can undergo interpolation
  std::unique_ptr<oops::FieldSet4D> fset4dXb;
  std::unique_ptr<oops::FieldSet4D> fset4dFg;

  // Change resolution if needed
  if (params.changeBackgroundResolution) {
    const State4D_ xb_lowres(geom, xb);
    const State4D_ fg_lowres(geom, fg);
    const oops::FieldSet4D fset4dXbTmp(xb_lowres);
    const oops::FieldSet4D fset4dFgTmp(fg_lowres);
    fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
    fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
  } else {
    const oops::FieldSet4D fset4dXbTmp(xb);
    const oops::FieldSet4D fset4dFgTmp(fg);
    fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
    fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
  }

  // Initialize outer variables
  const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
  oops::Variables outerVars(incVars);
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    outerVars[i].setLevels(vlevs[i]);
  }

  // Create blockchain
  blockChain_ = SaberBlockChainFactory<MODEL>::create(
        geom,
        outerVars,
        *fset4dXb,
        *fset4dFg,
        config);

  oops::Log::trace() << "ErrorCovariance::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// ----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doRandomize(Increment4D_ & dx) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  // This extra fieldset is only needed for backward compatibility in tests:
  // this allows the fields in the final fieldset to be in the same order as
  // blockChain_->outerVariables()
  oops::FieldSet4D fset4dSum(dx.times(), dx.commTime(), dx.geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
    fset4dSum[jtime].init(blockChain_->outerFunctionSpace(),
                          blockChain_->outerVariables(),
                          0.0);
  }

  // Create FieldSet4D, run randomize on it
  oops::FieldSet4D fset4d(dx.times(), dx.commTime(), dx.geometry().getComm());

  // Draw a random sample from the covariance
  blockChain_->randomize(fset4d);

  // For backward compatibility in tests
  fset4dSum += fset4d;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dx.size(); ++jtime) {
    dx[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment4D_ & dxi,
                                        Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  dxo = dxi;
  oops::FieldSet4D fset4d(dxo);

  // Extra fieldset only needed for backward compatibility in tests
  // that allows the fields in the final fieldset to be in the same order as
  // blockChain_->outerVariables()
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4d);
  fset4dSum.zero();

  // Apply covariance multiplication
  blockChain_->multiply(fset4d);

  // For backward compatibility in tests
  fset4dSum += fset4d;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dxo.size(); ++jtime) {
    dxo[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  os << "SABER ErrorCovariance";
  oops::Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
