/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/StdDev.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/oops/Utilities.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<StdDev> makerStdDev_("BUMP_StdDev");

// -----------------------------------------------------------------------------

StdDev::StdDev(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

  // Read input ATLAS files
  bump_->readAtlasFiles();

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

StdDev::~StdDev() {
  oops::Log::trace() << classname() << "::~StdDev starting" << std::endl;
  util::Timer timer(classname(), "~StdDev");
  oops::Log::trace() << classname() << "::~StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyStdDev(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyStdDev(fset);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> StdDev::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;
  std::vector<eckit::LocalConfiguration> inputModelFilesConf
    = bumpParams_.inputModelFilesConf.value().get_value_or({});
  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return bump_->getReadConfs(inputModelFilesConf);
}

// -----------------------------------------------------------------------------

void StdDev::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;
  for (const auto & fset : fsetVec) {
    bump_->addField(fset);
  }
  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  bump_->addEnsemble(fsetEns);
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;
  memberIndex_ = 0;
  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;
  bump_->iterativeUpdate(fset, memberIndex_);
  ++memberIndex_;
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationFinal() {
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;
  bump_->writeAtlasFiles();
  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> StdDev::fieldsToWrite() const {
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

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
