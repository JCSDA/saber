/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensity.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensity> makerDryAirDensity_("mo_dry_air_density");

// -----------------------------------------------------------------------------

DryAirDensity::DryAirDensity(const oops::GeometryData & outputGeometryData,
                             const std::vector<size_t> & activeVariableSizes,
                             const oops::Variables & outputVars,
                             const Parameters_ & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const std::vector<atlas::FieldSet> & fsetVec)
  : inputGeometryData_(outputGeometryData), inputVars_(outputVars), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::DryAirDensity starting" << std::endl;

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredStateVariables{ "exner_levels_minus_one",
                                                   "potential_temperature",
                                                   "exner",
                                                   "air_pressure_levels_minus_one",
                                                   "air_temperature",
                                                   "dry_air_density_levels_minus_one"};

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  for (auto & s : requiredStateVariables) {
    if (!xb.has(s)) {
      oops::Log::error() << "::DryAirDensity variable " << s <<
                            "is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(outputGeometryData.fieldSet()[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalDryAirDensity(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensity::~DryAirDensity() {
  oops::Log::trace() << classname() << "::~DryAirDensity starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensity");
  oops::Log::trace() << classname() << "::~DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalDryAirDensityTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------

void DryAirDensity::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalDryAirDensityAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
