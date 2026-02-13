/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberBlockParametersBase.h"

#include "eckit/exception/Exceptions.h"

namespace saber {

// -----------------------------------------------------------------------------

bool SaberBlockParametersBase::doCalibration() const {
  const auto & readConf = readParams.value();
  const auto & calibrationConf = calibrationParams.value();
  if (readConf != boost::none && calibrationConf != boost::none) {
    throw eckit::UserError("read and calibration configuration cannot be present at the same time"
      " in block " + saberBlockName.value(), Here());
  }
  const auto & blockName = saberBlockName.value();
  return (calibrationConf != boost::none || blockName == "Ensemble");
}

// -----------------------------------------------------------------------------

bool SaberBlockParametersBase::doRead() const {
  const auto & readConf = readParams.value();
  const auto & calibrationConf = calibrationParams.value();
  if (readConf != boost::none && calibrationConf != boost::none) {
    throw eckit::UserError("read and calibration configuration cannot be present at the same time"
      " in block " + saberBlockName.value(), Here());
  }
  return (readConf != boost::none);
}

// -----------------------------------------------------------------------------

oops::Variables SaberBlockParametersBase::getActiveVars(const oops::Variables & defaultVars) const {
  oops::Variables activeVars_nomd;
  if (this->mandatoryActiveVars().size() == 0) {
    // No mandatory active variables for this block
    activeVars_nomd = this->activeVars.value().get_value_or(defaultVars);
  } else {
    // Block with mandatory active variables
    activeVars_nomd = this->activeVars.value().get_value_or(this->mandatoryActiveVars());
    ASSERT(this->mandatoryActiveVars() <= activeVars_nomd);
  }
  // Copy the variables that exist in defaultVars from defaultVars (they have metadata
  // associated with them)
  oops::Variables activeVars;
  for (auto & var : activeVars_nomd) {
    if (defaultVars.has(var.name())) {
      activeVars.push_back(defaultVars[var.name()]);
    } else {
      activeVars.push_back(var);
    }
  }
  return activeVars;
}

// -----------------------------------------------------------------------------

}  // namespace saber
