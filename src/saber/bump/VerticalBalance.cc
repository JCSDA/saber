/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/VerticalBalance.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/oops/Utilities.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VerticalBalance> makerVerticalBalance_("BUMP_VerticalBalance");

// -----------------------------------------------------------------------------

VerticalBalance::VerticalBalance(const oops::GeometryData & outerGeometryData,
                                 const oops::Variables & outerVars,
                                 const eckit::Configuration & covarConf,
                                 const Parameters_ & params,
                                 const oops::FieldSet3D & xb,
                                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    bumpParams_(),
    bump_(),
    memberIndex_(0)
{
  oops::Log::trace() << classname() << "::VerticalBalance starting"
                     << std::endl;

  // Get active variables
  activeVars_ = getActiveVars(params, outerVars);
  std::vector<size_t> activeVariableSizes;
  for (std::string var : activeVars_.variables()) {
    activeVariableSizes.push_back(activeVars_.getLevels(var));
  }

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
                                 activeVars_.variables(),
                                 covarConf,
                                 bumpParams_.toConfiguration()));

  oops::Log::trace() << classname() << "::VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

VerticalBalance::~VerticalBalance() {
  oops::Log::trace() << classname() << "::~VerticalBalance starting" << std::endl;
  util::Timer timer(classname(), "~VerticalBalance");
  oops::Log::trace() << classname() << "::~VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyVbal(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyVbal(fset);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::directCalibration(const std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  bump_->addEnsemble(fsetEns);
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;
  memberIndex_ = 0;
  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::iterativeCalibrationUpdate(const atlas::FieldSet & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;
  bump_->iterativeUpdate(fset, memberIndex_);
  ++memberIndex_;
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::iterativeCalibrationFinal() {
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;
  bump_->writeAtlasFiles();
  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> VerticalBalance::fieldsToWrite()
  const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = bumpParams_.outputModelFilesConf.value().get_value_or({});
  if (outputModelFilesConf.size() > 0) {
    // Print log info
    oops::Log::info() << "Info     : "
                      << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      << std::endl;
    oops::Log::info() << "Info     : +++ Get parameters from BUMP"
                      << std::endl;
  }

  // Return configuration/fieldset pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> pairs
  = bump_->fieldsToWrite(outputModelFilesConf);

  if (outputModelFilesConf.size() > 0) {
    // Print log info
    oops::Log::info() << "Info     : "
                      << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      << std::endl;
    oops::Log::info() << "Info     : +++ Write files"
                      << std::endl;
  }
  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void VerticalBalance::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
