/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_BUMP_PSICHITOUV_H_
#define SABER_OOPS_BUMP_PSICHITOUV_H_

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
class BUMP_PsiChiToUVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_PsiChiToUVParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters<MODEL>> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class BUMP_PsiChiToUV : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;
  typedef BUMP<MODEL>           BUMP_;

 public:
  static const std::string classname() {return "saber::BUMP_PsiChiToUV";}

  typedef BUMP_PsiChiToUVParameters<MODEL> Parameters_;

  BUMP_PsiChiToUV(const Geometry_ &,
                  const Parameters_ & params,
                  const State_ &,
                  const State_ &);
  virtual ~BUMP_PsiChiToUV();

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
BUMP_PsiChiToUV<MODEL>::BUMP_PsiChiToUV(const Geometry_ & resol,
                                        const Parameters_ & params,
                                        const State_ & xb,
                                        const State_ & fg)
  : SaberBlockBase<MODEL>(params), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_PsiChiToUV starting" << std::endl;

  // Setup and check input/ouput variables (only two variables should be different)
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  size_t sameVarsInput = 0;
  for (size_t jvar = 0; jvar < inputVars.size(); ++jvar) {
    if (outputVars.has(inputVars[jvar])) sameVarsInput += 1;
  }
  ASSERT(sameVarsInput == inputVars.size()-2);
  size_t sameVarsOutput = 0;
  for (size_t jvar = 0; jvar < outputVars.size(); ++jvar) {
    if (inputVars.has(outputVars[jvar])) sameVarsOutput += 1;
  }
  ASSERT(sameVarsOutput == outputVars.size()-2);

  // Active variables
  oops::Variables allVars;
  allVars += inputVars;
  allVars += outputVars;
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= allVars);
  } else {
    activeVars += allVars;
  }

  // Initialize BUMP
  bump_.reset(new BUMP_(resol, activeVars, params.bumpParams.value(), xb, fg));

  oops::Log::trace() << classname() << "::BUMP_PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
BUMP_PsiChiToUV<MODEL>::~BUMP_PsiChiToUV() {
  oops::Log::trace() << classname() << "::~BUMP_PsiChiToUV starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_PsiChiToUV");
  oops::Log::trace() << classname() << "::~BUMP_PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(atlasFieldSet);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyPsiChiToUV(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet)
  const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  ABORT("BUMP_PsiChiToUV<MODEL>::inverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet)
  const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("BUMP_PsiChiToUV<MODEL>::inverseMultiplyAD: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void BUMP_PsiChiToUV<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_BUMP_PSICHITOUV_H_
