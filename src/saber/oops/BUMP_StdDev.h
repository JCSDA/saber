/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_BUMP_STDDEV_H_
#define SABER_OOPS_BUMP_STDDEV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"

#include "saber/oops/BUMP.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class BUMP_StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_StdDevParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------


template <typename MODEL>
class BUMP_StdDev : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef BUMP<MODEL>           BUMP_;

 public:
  static const std::string classname() {return "saber::BUMP_StdDev";}

  typedef BUMP_StdDevParameters Parameters_;

  BUMP_StdDev(const Geometry_ &, const Parameters_ &);
  virtual ~BUMP_StdDev();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<BUMP_> bump_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_StdDev<MODEL>::BUMP_StdDev(const Geometry_ & resol, const Parameters_ & params)
  : SaberBlockBase<MODEL>(params), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_StdDev starting" << std::endl;

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
  bump_.reset(new BUMP_(resol, activeVars, params.bumpParams.value()));

  oops::Log::trace() << classname() << "::BUMP_StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_StdDev<MODEL>::~BUMP_StdDev() {
  oops::Log::trace() << classname() << "::~BUMP_StdDev starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_StdDev");
  oops::Log::trace() << classname() << "::~BUMP_StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(atlasFieldSet);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyStdDev(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  bump_->inverseMultiplyStdDev(atlasFieldSet);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  this->inverseMultiply(atlasFieldSet);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_StdDev<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_BUMP_STDDEV_H_
