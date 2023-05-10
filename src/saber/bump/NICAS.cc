/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/NICAS.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/bump/lib/Utilities.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<NICAS> makerNICAS_("BUMP_NICAS");

// -----------------------------------------------------------------------------

NICAS::NICAS(const oops::GeometryData & geometryData,
             const std::vector<size_t> & activeVariableSizes,
             const oops::Variables & centralVars,
             const eckit::Configuration & covarConf,
             const Parameters_ & params,
             const atlas::FieldSet & xb,
             const atlas::FieldSet & fg,
             const util::DateTime & validTimeOfXbFg,
             const size_t & timeRank)
  : SaberCentralBlockBase(params), bumpParams_(),
    bump_(),
    memberIndex_(0)
{
  oops::Log::trace() << classname() << "::NICAS starting" << std::endl;

  // Get active variables
  activeVars_ = params.activeVars.value().get_value_or(centralVars);

  // Get BUMP parameters
  if (params.doCalibration()) {
    bumpParams_ = *params.calibrationParams.value();
  } else if (params.doRead()) {
    bumpParams_ = *params.readParams.value();
  } else {
    ABORT("calibration or read required in BUMP");
  }

  // Initialize BUMP
  bump_.reset(new bump_lib::BUMP(geometryData.comm(),
                                 oops::LibOOPS::instance().infoChannel(),
                                 oops::LibOOPS::instance().testChannel(),
                                 geometryData.functionSpace(),
                                 geometryData.fieldSet(),
                                 activeVariableSizes,
                                 activeVars_.variables(),
                                 covarConf,
                                 bumpParams_.toConfiguration(),
                                 timeRank));

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

void NICAS::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  bump_->randomizeNicas(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyNicas(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> NICAS::fieldsToRead() {
  oops::Log::trace() << classname() << "::fieldsToRead starting" << std::endl;
  std::vector<eckit::LocalConfiguration> inputModelFilesConf
    = bumpParams_.inputModelFilesConf.value().get_value_or({});
  return bump_->fieldsToRead(inputModelFilesConf);
}

// -----------------------------------------------------------------------------

void NICAS::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  for (const auto & input : bump_->inputs()) {
    bump_->addField(input.second);
  }
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::directCalibration(const std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  for (const auto & input : bump_->inputs()) {
    bump_->addField(input.second);
  }
  bump_->addEnsemble(fsetEns);
  bump_->runDrivers();
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;
  for (const auto & input : bump_->inputs()) {
    bump_->addField(input.second);
  }
  memberIndex_ = 0;
  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::iterativeCalibrationUpdate(const atlas::FieldSet & fset) {
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

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> NICAS::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = bumpParams_.outputModelFilesConf.value().get_value_or({});
  return bump_->fieldsToWrite(outputModelFilesConf);
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
