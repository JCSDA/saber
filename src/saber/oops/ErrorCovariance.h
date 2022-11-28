/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
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
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralTBlock.h"
#include "saber/oops/SaberOuterTBlock.h"

namespace saber {

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<MODEL>)
 public:
  oops::RequiredParameter<SaberCentralTBlockParameters<MODEL>>
    saberCentralTBlockParams{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterTBlockParameters<MODEL>>>
    saberOuterTBlocksParams{"saber outer blocks", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                 Geometry_;
  typedef oops::Increment<MODEL>                                Increment_;
  typedef typename boost::ptr_vector<SaberOuterTBlock<MODEL>>   SaberOuterTBlockVec_;
  typedef typename SaberOuterTBlockVec_::iterator               iter_;
  typedef typename SaberOuterTBlockVec_::const_iterator         icst_;
  typedef typename SaberOuterTBlockVec_::const_reverse_iterator ircst_;
  typedef oops::State<MODEL>                                    State_;
  typedef typename std::unordered_map<std::string, const oops::GeometryData*> GeometryDataMap_;

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

  GeometryDataMap_ outerGeometryDataMap_;
  std::unique_ptr<SaberCentralTBlock<MODEL>> saberCentralTBlock_;
  SaberOuterTBlockVec_ saberOuterTBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const Parameters_ & params,
                                        const State_ & xb,
                                        const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, params, xb, fg), saberCentralTBlock_(),
    saberOuterTBlocks_()
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Local copy of background and first guess
  State_ xbLocal(xb);
  State_ fgLocal(fg);

  // Initialize geometryData map
  for (const auto var : incVars.variables()) {
    outerGeometryDataMap_[var] = &(geom.generic());
  }

  // Intialize outer variables
  oops::Variables outerVars(incVars);

  // Build outer blocks successively
  const boost::optional<std::vector<SaberOuterTBlockParameters<MODEL>>> &saberOuterTBlocksParams =
    params.saberOuterTBlocksParams.value();
  if (saberOuterTBlocksParams != boost::none) {
    // Loop in reverse order
    for (const SaberOuterTBlockParameters<MODEL> & saberOuterTBlockParams :
      boost::adaptors::reverse(*saberOuterTBlocksParams)) {
      // Create outer templated block
      saberOuterTBlocks_.push_back(new SaberOuterTBlock<MODEL>(geom,
                                   outerGeometryDataMap_,
                                   outerVars,
                                   saberOuterTBlockParams,
                                   xbLocal,
                                   fgLocal));

      // Update outer geometry and variables for the next block
      outerVars = saberOuterTBlocks_.back().innerVars();
      outerGeometryDataMap_ = saberOuterTBlocks_.back().innerGeometryDataMap();
    }
  }

  // Get central templated block parameters
  const SaberCentralTBlockParameters<MODEL> saberCentralTBlockParams =
    params.saberCentralTBlockParams.value();

  // Create central templated block
  saberCentralTBlock_.reset(new SaberCentralTBlock<MODEL>(geom,
                            outerGeometryDataMap_,
                            outerVars,
                            saberCentralTBlockParams,
                            xbLocal,
                            fgLocal));

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

  // Random vector (necessary for some SABER blocks)
  dx.random();

  // Central block randomization
  saberCentralTBlock_->randomize(dx.fieldSet());

  // Outer blocks forward multiplication
  for (ircst_ it = saberOuterTBlocks_.rbegin(); it != saberOuterTBlocks_.rend(); ++it) {
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
  for (icst_ it = saberOuterTBlocks_.begin(); it != saberOuterTBlocks_.end(); ++it) {
    it->multiplyAD(dxo.fieldSet());
  }

  // Central block multiplication
  if (saberCentralTBlock_) {
    saberCentralTBlock_->multiply(dxo.fieldSet());
  }

  // Outer blocks forward multiplication
  for (ircst_ it = saberOuterTBlocks_.rbegin(); it != saberOuterTBlocks_.rend(); ++it) {
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
