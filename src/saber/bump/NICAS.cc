/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/NICAS.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/oops/Utilities.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<NICAS> makerNICAS_("BUMP_NICAS");

// -----------------------------------------------------------------------------

NICAS::NICAS(const oops::GeometryData & geometryData,
             const oops::Variables & centralVars,
             const eckit::Configuration & covarConf,
             const Parameters_ & params,
             const oops::FieldSet3D & xb,
             const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()),
    activeVars_(getActiveVars(params, centralVars)),
    bumpParams_(params.calibrationParams.value() != boost::none ? *params.calibrationParams.value()
      : *params.readParams.value()),
    bump_(new BUMP(geometryData, activeVars_, covarConf, bumpParams_,
      params.fieldsMetaData.value(), xb)),
    memberIndex_(0) {
  oops::Log::trace() << classname() << "::NICAS starting" << std::endl;

  // Read input ATLAS files
  bump_->readAtlasFiles();

  oops::Log::trace() << classname() << "::NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

NICAS::~NICAS() {
  oops::Log::trace() << classname() << "::~NICAS starting" << std::endl;
  util::Timer timer(classname(), "~NICAS");
  oops::Log::trace() << classname() << "::~NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  bump_->randomizeNicas(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyNicas(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> NICAS::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;
  std::vector<eckit::LocalConfiguration> inputModelFilesConf
    = bumpParams_.inputModelFilesConf.value().get_value_or({});
  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return bump_->getReadConfs(inputModelFilesConf);
}

// -----------------------------------------------------------------------------

void NICAS::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;
  for (const auto & fset : fsetVec) {
    bump_->addField(fset);
  }
  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  bump_->addEnsemble(fsetEns);
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;
  memberIndex_ = 0;
  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;
  bump_->iterativeUpdate(fset, memberIndex_);
  ++memberIndex_;
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::iterativeCalibrationFinal() {
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal starting" << std::endl;
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::dualResolutionSetup(const oops::GeometryData & geometryData) {
  oops::Log::trace() << classname() << "::dualResolutionSetup starting" << std::endl;
  bump_->dualResolutionSetup(geometryData.functionSpace(),
                             geometryData.fieldSet());
  oops::Log::trace() << classname() << "::dualResolutionSetup done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::multiplySqrt(const atlas::Field & cv,
                         oops::FieldSet3D & fset,
                         const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;
  bump_->multiplyNicasSqrt(cv, fset, offset);
  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::multiplySqrtAD(const oops::FieldSet3D & fset,
                           atlas::Field & cv,
                           const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;
  bump_->multiplyNicasSqrtAd(fset, cv, offset);
  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> NICAS::fieldsToWrite() const {
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

void NICAS::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;
  bump_->writeAtlasFiles();
  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
