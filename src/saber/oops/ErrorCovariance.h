/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
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
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

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
    saberCentralBlock{"saber central block", this};
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

  std::vector<std::reference_wrapper<const oops::GeometryData>> outerGeometryData_;
  std::unique_ptr<SaberCentralBlockBase> saberCentralBlock_;
  SaberOuterBlockVec_ saberOuterBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const Parameters_ & params,
                                        const State_ & xb,
                                        const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, params, xb, fg), saberCentralBlock_(),
    saberOuterBlocks_()
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Local copy of background and first guess
  State_ xbLocal(xb);
  State_ fgLocal(fg);

  // Initial outer geometry and variables
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
      std::vector<atlas::FieldSet> fsetVec = readInputFields(
        geom,
        activeVars,
        xb.validTime(),
        saberOuterBlockParams.inputFields.value());

      // Create outer block
      oops::Log::info() << "Creating outer block: " << saberOuterBlockParams.saberBlockName.value()
                        << std::endl;
      saberOuterBlocks_.push_back(SaberOuterBlockFactory::create(
                                  outerGeometryData_.back().get(),
                                  geom.variableSizes(activeVars),
                                  outerVars,
                                  saberOuterBlockParams,
                                  xbLocal.fieldSet(),
                                  fgLocal.fieldSet(),
                                  fsetVec));

      // Apply calibration inverse on xb and fg
      // TODO(Benjamin): uncomment these lines when all blocks are compliant
//      saberOuterBlocks_.back().calibrationInverseMultiply(xbLocal.fieldSet());
//      saberOuterBlocks_.back().calibrationInverseMultiply(fgLocal.fieldSet());

      // Access inner geometry and variables
      const oops::GeometryData & innerGeometryData = saberOuterBlocks_.back().innerGeometryData();
      const oops::Variables innerVars = saberOuterBlocks_.back().innerVars();

      // Check that active variables are present in either inner or outer variables, or both
      for (const auto & var : activeVars.variables()) {
        ASSERT(innerVars.has(var) || outerVars.has(var));
      }

      // Update outer geometry and variables for the next block
      outerGeometryData_.push_back(innerGeometryData);
      outerVars = innerVars;
    }
  }

  // Get central block parameters
  const SaberCentralBlockParametersWrapper & saberCentralBlockParamWrapper =
    params.saberCentralBlock.value();
  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamWrapper.saberCentralBlockParameters;

  // Define centra variables
  oops::Variables centralVars = outerVars;

  // Get active variables
  oops::Variables activeVars = saberCentralBlockParams.activeVars.value().get_value_or(outerVars);

  // Read input fields (on model increment geometry)
  std::vector<atlas::FieldSet> fsetVec = readInputFields(
    geom,
    activeVars,
    xb.validTime(),
    saberCentralBlockParams.inputFields.value());

  // Create central block
  saberCentralBlock_.reset(SaberCentralBlockFactory::create(
                           outerGeometryData_.back().get(),
                           geom.variableSizes(activeVars),
                           centralVars,
                           saberCentralBlockParams,
                           xbLocal.fieldSet(),
                           fgLocal.fieldSet(),
                           fsetVec));

  // Check that active variables are present in central variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(centralVars.has(var));
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

  // Random vector (necessary for some SABER blocks)
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
