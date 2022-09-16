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
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BUMP_PsiChiToUV> makerBUMP_PsiChiToUV_("BUMP_PsiChiToUV");

// -----------------------------------------------------------------------------

BUMP_PsiChiToUV::BUMP_PsiChiToUV(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_PsiChiToUV starting" << std::endl;

  // Deserialize configuration
  BUMP_PsiChiToUVParameters params;
  params.deserialize(conf);

  // Input geometry and variables
  inputFunctionSpace_ = outputFunctionSpace;
  inputExtraFields_ = outputExtraFields;
  inputVars_ = params.outputVars.value(); // TODO(Benjamin): remove u/v, add psi/chi

  // Initialize BUMP
  bump_.reset(new BUMP(comm,
                       outputFunctionSpace,
                       outputExtraFields,
                       activeVariableSizes,
                       *params.activeVars.value(),
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

void BUMP_PsiChiToUV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyPsiChiToUV(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::calibrationInverseMultiply(atlas::FieldSet & fset)
  const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_PsiChiToUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
