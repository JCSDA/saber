/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/StdDev.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<StdDev> makerStdDev_("StdDev");

// -----------------------------------------------------------------------------

StdDev::StdDev(const oops::GeometryData & outerGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const oops::Variables & outerVars,
               const eckit::Configuration & covarConfig,
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData),
    activeVariableSizes_(activeVariableSizes),
    innerVars_(outerVars),
    params_(params),
    readFromAtlas_(false),
    readFromModel_(false),
    writeToAtlas_(false),
    writeToModel_(false)
{
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

  // Get active variables
  oops::Variables activeVars = params_.activeVars.value().get_value_or(outerVars);

  // Prepare read parameters
  const auto & readParams = params_.readParams.value();
  if (readParams != boost::none) {
    const auto & atlasFileConf = readParams->atlasFileConf.value();
    readFromAtlas_ = (atlasFileConf != boost::none);
    const auto & modelFileConf = readParams->modelFileConf.value();
    readFromModel_ = (modelFileConf != boost::none);
    if (readFromAtlas_ && readFromModel_) {
      ABORT("either ATLAS or model stddev input file should be present, not both");
    }
    if (readFromAtlas_) {
      readConf_ = *atlasFileConf;
    }
    if (readFromModel_) {
      readConf_ = *modelFileConf;
    }
  }

  // Prepare write parameters
  const auto & calibrationParams = params_.calibrationParams.value();
  if (calibrationParams != boost::none) {
    const auto & atlasFileConf = calibrationParams->atlasFileConf.value();
    writeToAtlas_ = (atlasFileConf != boost::none);
    const auto & modelFileConf = calibrationParams->modelFileConf.value();
    writeToModel_ = (modelFileConf != boost::none);
    if (writeToAtlas_ && writeToModel_) {
      ABORT("either ATLAS or model stddev output file should be present, not both");
    }
    if (writeToAtlas_) {
      writeConf_ = *atlasFileConf;
    }
    if (writeToModel_) {
      writeConf_ = *modelFileConf;
    }
  }

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::multiplyFieldSets(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::multiplyFieldSets(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  util::divideFieldSets(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> StdDev::fieldsToRead() {
  oops::Log::trace() << classname() << "::fieldsToRead starting" << std::endl;

  if (readFromModel_) {
    // Create FieldSet
    atlas::FieldSet fset;
    fset.name() = "StdDev";

    // Add pair
    inputs_.push_back(std::make_pair(readConf_, fset));
  }

  oops::Log::trace() << classname() << "::fieldsToRead done" << std::endl;
  return inputs_;
}

// -----------------------------------------------------------------------------

void StdDev::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Read ATLAS stddev file
  if (readFromAtlas_) {
    // Read file
    util::readFieldSet(innerGeometryData_.comm(),
                       innerGeometryData_.functionSpace(),
                       activeVariableSizes_,
                       innerVars_.variables(),
                       readConf_,
                       stdDevFset_);

    // Set name
    stdDevFset_.name() = "StdDev";

    // Print FieldSet norm
    oops::Log::test() << "Norm of input parameter StdDev: "
                      << util::normFieldSet(stdDevFset_,
                                            innerVars_.variables(),
                                            innerGeometryData_.comm())
                      << std::endl;
  }

  // Use model stddev file
  if (readFromModel_) {
    // Copy file
    util::copyFieldSet(inputs_[0].second, stdDevFset_);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::directCalibration(const std::vector<atlas::FieldSet> & fsetEns) {
  // Initialize
  atlas::FieldSet mean;
  atlas::FieldSet var;
  for (size_t jvar = 0; jvar < innerVars_.size(); ++jvar) {
    mean.add(innerGeometryData_.functionSpace().createField<double>(
             atlas::option::name(innerVars_[jvar]) |
             atlas::option::levels(activeVariableSizes_[jvar])));
    var.add(innerGeometryData_.functionSpace().createField<double>(
            atlas::option::name(innerVars_[jvar]) |
            atlas::option::levels(activeVariableSizes_[jvar])));
  }
  util::zeroFieldSet(mean);
  util::zeroFieldSet(var);

  // Compute mean
  for (const auto & fset : fsetEns) {
    util::addFieldSets(mean, fset);
  }
  const double facMean = 1.0/static_cast<double>(fsetEns.size());
  util::multiplyFieldSet(mean, facMean);

  // Compute variance
  for (const auto & fset : fsetEns) {
    // Compute perturbation
    atlas::FieldSet fset_pert = util::copyFieldSet(fset);
    util::subtractFieldSets(fset_pert, mean);

    // Compute squared pertrubation
    atlas::FieldSet fset_pert_squared = util::copyFieldSet(fset_pert);
    util::multiplyFieldSets(fset_pert_squared, fset_pert);

    // Update sum
    util::addFieldSets(var, fset_pert_squared);
  }
  const double facVar = 1.0/static_cast<double>(fsetEns.size()-1);
  util::multiplyFieldSet(var, facVar);

  // Get standard deviation
  stdDevFset_ = util::copyFieldSet(var);
  util::sqrtFieldSet(stdDevFset_);

  // Set name
  stdDevFset_.name() = "StdDev";
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationInit() {
  // Initialize iterative counters with zeroes
  iterativeN_ = 0;
  iterativeMean_.clear();
  iterativeVar_.clear();
  for (size_t jvar = 0; jvar < innerVars_.size(); ++jvar) {
    iterativeMean_.add(innerGeometryData_.functionSpace().createField<double>(
                       atlas::option::name(innerVars_[jvar]) |
                       atlas::option::levels(activeVariableSizes_[jvar])));
    iterativeVar_.add(innerGeometryData_.functionSpace().createField<double>(
                      atlas::option::name(innerVars_[jvar]) |
                      atlas::option::levels(activeVariableSizes_[jvar])));
  }
  util::zeroFieldSet(iterativeMean_);
  util::zeroFieldSet(iterativeVar_);
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationUpdate(const atlas::FieldSet & fset) {
  // Increment ensemble index
  // ie = ie + 1
  iterativeN_++;

  // Remove mean
  // pert = state - mean
  atlas::FieldSet fset_pert = util::copyFieldSet(fset);
  util::subtractFieldSets(fset_pert, iterativeMean_);

  // Update variance
  // var = var + (ie-1)/ie * pert^2
  atlas::FieldSet fset_pert_squared = util::copyFieldSet(fset_pert);
  util::multiplyFieldSets(fset_pert_squared, fset_pert);
  const double facVar = static_cast<double>(iterativeN_-1)/static_cast<double>(iterativeN_);
  util::multiplyFieldSet(fset_pert_squared, facVar);
  util::addFieldSets(iterativeVar_, fset_pert_squared);

  // Update mean
  // mean = mean + 1 / ie * pert
  const double facMean = 1.0/static_cast<double>(iterativeN_);
  util::multiplyFieldSet(fset_pert, facMean);
  util::addFieldSets(iterativeMean_, fset_pert);
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationFinal() {
  // Normalize variance
  // var = 1 / (N-1) * var
  const double facVar = 1.0/static_cast<double>(iterativeN_-1);
  util::multiplyFieldSet(iterativeVar_, facVar);

  // Get standard-deviation
  stdDevFset_ = util::copyFieldSet(iterativeVar_);
  util::sqrtFieldSet(stdDevFset_);

  // Set name
  stdDevFset_.name() = "StdDev";

  // Cleaning
  iterativeN_ = 0;
  iterativeMean_.clear();
  iterativeVar_.clear();
}

// -----------------------------------------------------------------------------

void StdDev::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Write ATLAS stddev file
  if (writeToAtlas_) {
    // Print FieldSet norm
    oops::Log::test() << "Norm of output parameter StdDev: "
                      << util::normFieldSet(stdDevFset_,
                                            innerVars_.variables(),
                                            innerGeometryData_.comm())
                      << std::endl;


    // Write file
    util::writeFieldSet(innerGeometryData_.comm(), writeConf_, stdDevFset_);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> StdDev::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> outputs;

  if (writeToModel_) {
    // Add pair
    outputs.push_back(std::make_pair(writeConf_, stdDevFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
