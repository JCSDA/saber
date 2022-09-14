/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP_PsiChiToUV.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/bump/BUMP.h"
#include "saber/bump/BUMP_Parameters.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<BUMP_PsiChiToUV> makerBUMP_PsiChiToUV_("BUMP_PsiChiToUV");

// -----------------------------------------------------------------------------

BUMP_PsiChiToUV::BUMP_PsiChiToUV(const eckit::mpi::Comm & comm,
                                 const atlas::FunctionSpace & functionSpace,
                                 const atlas::FieldSet & extraFields,
                                 const std::vector<size_t> & variableSizes,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), bump_()
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
  bump_.reset(new BUMP(comm,
                       functionSpace,
                       extraFields,
                       variableSizes,
                       activeVars,
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::BUMP_PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP_PsiChiToUV::~BUMP_PsiChiToUV() {
  oops::Log::trace() << classname() << "::~BUMP_PsiChiToUV starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_PsiChiToUV");
  oops::Log::trace() << classname() << "::~BUMP_PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyPsiChiToUV(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::inverseMultiply(atlas::FieldSet & fset)
  const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  ABORT("BUMP_PsiChiToUV::inverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::inverseMultiplyAD(atlas::FieldSet & fset)
  const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("BUMP_PsiChiToUV::inverseMultiplyAD: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
