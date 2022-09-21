/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensitySaberBlock.h"

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

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensitySaberBlock>
       makerDryAirDensitySaberBlock_("mo_dry_air_density");

// -----------------------------------------------------------------------------

DryAirDensitySaberBlock::DryAirDensitySaberBlock(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::DryAirDensitySaberBlock starting" << std::endl;

  // Deserialize configuration
  DryAirDensitySaberBlockParameters params;
  params.deserialize(conf);

  // Input geometry and variables
  inputFunctionSpace_ = outputFunctionSpace;
  inputExtraFields_ = outputExtraFields;
  inputVars_ = params.outputVars.value();

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
    if (!xb.has_field(s)) {
      oops::Log::error() << "::DryAirDensitySaberBlock variable " << s <<
                            "is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(outputExtraFields[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalDryAirDensity(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::DryAirDensitySaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensitySaberBlock::~DryAirDensitySaberBlock() {
  oops::Log::trace() << classname() << "::~DryAirDensitySaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensitySaberBlock");
  oops::Log::trace() << classname() << "::~DryAirDensitySaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalDryAirDensityTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalDryAirDensityAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------


}  // namespace saber
