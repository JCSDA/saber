/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_BUMP_VERTICALBALANCE_H_
#define SABER_BUMP_BUMP_VERTICALBALANCE_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
class BUMP_VerticalBalanceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_VerticalBalanceParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters<MODEL>> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class BUMP_VerticalBalance : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;
  typedef BUMP<MODEL>           BUMP_;

 public:
  static const std::string classname() {return "saber::BUMP_VerticalBalance";}

  typedef BUMP_VerticalBalanceParameters<MODEL> Parameters_;

  BUMP_VerticalBalance(const Geometry_ &,
                       const Parameters_ & params,
                       const State_ &,
                       const State_ &);
  virtual ~BUMP_VerticalBalance();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<BUMP_> bump_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_VerticalBalance<MODEL>::BUMP_VerticalBalance(const Geometry_ & resol,
                                                  const Parameters_ & params,
                                                  const State_ & xb,
                                                  const State_ & fg)
  : SaberBlockBase<MODEL>(params), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_VerticalBalance starting"
                     << std::endl;

  // Setup and check input/ouput variables
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  ASSERT(inputVars == outputVars);

  // Active variables
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= inputVars);
  } else {
    activeVars += inputVars;
  }

  // Initialize BUMP
  bump_.reset(new BUMP_(resol, activeVars, params.bumpParams.value(), xb, fg));

  oops::Log::trace() << classname() << "::BUMP_VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_VerticalBalance<MODEL>::~BUMP_VerticalBalance() {
  oops::Log::trace() << classname() << "::~BUMP_VerticalBalance starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_VerticalBalance");
  oops::Log::trace() << classname() << "::~BUMP_VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyVbal(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  bump_->inverseMultiplyVbal(fset);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  bump_->inverseMultiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_VerticalBalance<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_BUMP_BUMP_VERTICALBALANCE_H_
