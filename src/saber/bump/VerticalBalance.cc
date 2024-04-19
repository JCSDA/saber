/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/VerticalBalance.h"

#include "eckit/exception/Exceptions.h"

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
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    activeVars_(getActiveVars(params, outerVars)),
    bumpParams_(params.calibrationParams.value() != boost::none ? *params.calibrationParams.value()
      : *params.readParams.value()),
    bump_(new BUMP(outerGeometryData, activeVars_, covarConf, bumpParams_,
      params.fieldsMetaData.value(), xb)),
    memberIndex_(0) {
  oops::Log::trace() << classname() << "::VerticalBalance starting" << std::endl;
  oops::Log::trace() << classname() << "::VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

VerticalBalance::~VerticalBalance() {
  oops::Log::trace() << classname() << "::~VerticalBalance starting" << std::endl;
  util::Timer timer(classname(), "~VerticalBalance");
  oops::Log::trace() << classname() << "::~VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyVbal(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::leftInverseMultiply(oops::FieldSet3D & fset) const {
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

void VerticalBalance::directCalibration(const oops::FieldSets & fsetEns) {
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

void VerticalBalance::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
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

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> VerticalBalance::fieldsToWrite()
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
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> pairs
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
