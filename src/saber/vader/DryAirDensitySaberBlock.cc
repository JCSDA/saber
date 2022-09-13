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

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<DryAirDensitySaberBlock>
       makerDryAirDensitySaberBlock_("mo_dry_air_density");

// -----------------------------------------------------------------------------

DryAirDensitySaberBlock::DryAirDensitySaberBlock(const eckit::mpi::Comm & comm,
                                                 const atlas::FunctionSpace & functionSpace,
                                                 const atlas::FieldSet & extraFields,
                                                 const std::vector<size_t> & variableSizes,
                                                 const Parameters_ & params,
                                                 const atlas::FieldSet & xb,
                                                 const atlas::FieldSet & fg,
                                                 const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::DryAirDensitySaberBlock starting" << std::endl;

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
    augmentedStateFieldSet_.add(extraFields[s]);
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

void DryAirDensitySaberBlock::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("DryAirDensitySaberBlock::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalDryAirDensityTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalDryAirDensityAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiplyAD not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensitySaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------


}  // namespace saber
