/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/AirTemperatureSaberBlock.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/control2analysis_linearvarchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<AirTemperatureSaberBlock>
       makerAirTemperatureSaberBlock_("mo_air_temperature");

// -----------------------------------------------------------------------------

AirTemperatureSaberBlock::AirTemperatureSaberBlock(const eckit::mpi::Comm & comm,
               const oops::GeometryData & outputGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), inputGeometryData_(outputGeometryData), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::AirTemperatureSaberBlock starting" << std::endl;

  // Deserialize configuration
  AirTemperatureSaberBlockParameters params;
  params.deserialize(conf);

  // Input variables
  inputVars_ = params.outputVars.value();

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredStateVariables{"exner_levels_minus_one",
                                                  "potential_temperature"};

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has_field(s)) {
      oops::Log::info() << "AirTemperatureSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(outputGeometryData.fieldSet()[s]);
  }

  oops::Log::trace() << classname() << "::AirTemperatureSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

AirTemperatureSaberBlock::~AirTemperatureSaberBlock() {
  oops::Log::trace() << classname() << "::~AirTemperatureSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~AirTemperatureSaberBlock");
  oops::Log::trace() << classname() << "::~AirTemperatureSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalAirTemperatureTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalAirTemperatureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
