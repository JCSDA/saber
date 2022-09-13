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

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<AirTemperatureSaberBlock>
       makerAirTemperatureSaberBlock_("mo_air_temperature");

// -----------------------------------------------------------------------------

AirTemperatureSaberBlock::AirTemperatureSaberBlock(const atlas::FunctionSpace & functionSpace,
                                                   const atlas::FieldSet & extraFields,
                                                   const std::vector<size_t> & variableSizes,
                                                   const Parameters_ & params,
                                                   const atlas::FieldSet & xb,
                                                   const atlas::FieldSet & fg,
                                                   const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::AirTemperatureSaberBlock starting" << std::endl;

  // Setup and check input/output variables
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  ASSERT(inputVars == outputVars);

  // Active variables
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= inputVars);
  } else {
    activeVars += inputVars;
  }

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
    augmentedStateFieldSet_.add(extraFields[s]);
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

void AirTemperatureSaberBlock::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("AirTemperatureSaberBlock::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalAirTemperatureTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalAirTemperatureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiplyAD not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperatureSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
