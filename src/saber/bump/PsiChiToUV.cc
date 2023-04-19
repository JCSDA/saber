/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/PsiChiToUV.h"

#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<PsiChiToUV> makerPsiChiToUV_("BUMP_PsiChiToUV");

// -----------------------------------------------------------------------------

PsiChiToUV::PsiChiToUV(const oops::GeometryData & outerGeometryData,
                       const std::vector<size_t> & activeVariableSizes,
                       const oops::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const atlas::FieldSet & xb,
                       const atlas::FieldSet & fg)
  : innerGeometryData_(outerGeometryData),
    innerVars_(oops::Variables({"stream_function", "velocity_potential"})),
    outerVars_(outerVars),
    levels_(activeVariableSizes[0]),
    bumpParams_(),
    bump_(),
    memberIndex_(0)
{
  oops::Log::trace() << classname() << "::PsiChiToUV starting" << std::endl;

  // Get BUMP parameters
  if (params.doCalibration()) {
    bumpParams_ = *params.calibrationParams.value();
  } else if (params.doRead()) {
    bumpParams_ = *params.readParams.value();
  } else {
    ABORT("calibration or read required in BUMP");
  }

  // Initialize BUMP
  bump_.reset(new bump_lib::BUMP(outerGeometryData.comm(),
                                 oops::LibOOPS::instance().infoChannel(),
                                 oops::LibOOPS::instance().testChannel(),
                                 outerGeometryData.functionSpace(),
                                 outerGeometryData.fieldSet(),
                                 activeVariableSizes,
                                 params.mandatoryActiveVars().variables(),
                                 covarConf,
                                 bumpParams_.toConfiguration()));

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
  util::removeFieldsFromFieldSet(fset, innerVars_.variables());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyPsiChiToUVAd(fset);
  util::removeFieldsFromFieldSet(fset, outerVars_.variables());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  ABORT(classname() + "::leftInverseMultiply not implemented");
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void PsiChiToUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
