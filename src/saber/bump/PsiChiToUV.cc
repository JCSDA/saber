/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/PsiChiToUV.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<PsiChiToUV> makerPsiChiToUV_("BUMP_PsiChiToUV");

// -----------------------------------------------------------------------------

PsiChiToUV::PsiChiToUV(const oops::GeometryData & outputGeometryData,
                       const std::vector<size_t> & activeVariableSizes,
                       const oops::Variables & outputVars,
                       const Parameters_ & params,
                       const atlas::FieldSet & xb,
                       const atlas::FieldSet & fg,
                       const std::vector<atlas::FieldSet> & fsetVec)
  : inputGeometryData_(outputGeometryData), bump_()
{
  oops::Log::trace() << classname() << "::PsiChiToUV starting" << std::endl;

  // Check that active variables are present in parameters
  ASSERT(params.activeVars.value() != boost::none);

  // Get active variables
  oops::Variables activeVars = *params.activeVars.value();

  // Check active variables size
  ASSERT(activeVars.size() == 4);

  // Only two active variables should be part of output variables, other two are input variables
  size_t activeVarsInOutput = 0;
  for (const auto var : outputVars.variables()) {
    if (activeVars.has(var)) {
      activeVarsInOutput += 1;
    } else {
      inputVars_.push_back(var);
    }
  }
  ASSERT(activeVarsInOutput == 2);
  for (const auto var : activeVars.variables()) {
    if (!outputVars.has(var)) {
      inputVars_.push_back(var);
    }
  }

  // Initialize BUMP
  bump_.reset(new BUMP(outputGeometryData.comm(),
                       outputGeometryData.functionSpace(),
                       outputGeometryData.fieldSet(),
                       activeVariableSizes,
                       activeVars,
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

PsiChiToUV::~PsiChiToUV() {
  oops::Log::trace() << classname() << "::~PsiChiToUV starting" << std::endl;
  util::Timer timer(classname(), "~PsiChiToUV");
  oops::Log::trace() << classname() << "::~PsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyPsiChiToUV(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::calibrationInverseMultiply(atlas::FieldSet & fset)
  const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
