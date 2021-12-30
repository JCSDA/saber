/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_BUMP_NICAS_H_
#define SABER_OOPS_BUMP_NICAS_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/BUMP.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
class BUMP_NICASParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_NICASParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters<MODEL>> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class BUMP_NICAS : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef BUMP<MODEL>           BUMP_;

 public:
  static const std::string classname() {return "saber::BUMP_NICAS";}

  typedef BUMP_NICASParameters<MODEL> Parameters_;

  BUMP_NICAS(const Geometry_ &, const Parameters_ &);
  virtual ~BUMP_NICAS();

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
BUMP_NICAS<MODEL>::BUMP_NICAS(const Geometry_ & resol, const Parameters_ & params)
  : SaberBlockBase<MODEL>(params), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_NICAS starting" << std::endl;

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

  oops::Log::trace() << classname() << "::BUMP_NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_NICAS<MODEL>::~BUMP_NICAS() {
  oops::Log::trace() << classname() << "::~BUMP_NICAS starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_NICAS");
  oops::Log::trace() << classname() << "::~BUMP_NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  bump_->randomizeNicas(atlasFieldSet);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyNicas(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  ABORT("BUMP_NICAS<MODEL>::inverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  ABORT("BUMP_NICAS<MODEL>::multiplyAD: not implemented");
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("BUMP_NICAS<MODEL>::inverseMultiplyAD: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_NICAS<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_BUMP_NICAS_H_
