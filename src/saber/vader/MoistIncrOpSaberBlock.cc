/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/MoistIncrOpSaberBlock.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "../../vader/src/mo/common_varchange.h"
#include "../../vader/src/mo/control2analysis_linearvarchange.h"
#include "../../vader/src/mo/functions.h"
#include "../../vader/src/mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
static SaberBlockMaker<MoistIncrOpSaberBlock>
       makerMoistIncrOpSaberBlock_("mo_moistincrop");

// -----------------------------------------------------------------------------

MoistIncrOpSaberBlock::MoistIncrOpSaberBlock(const atlas::FunctionSpace & functionSpace,
                                             const atlas::FieldSet & extraFields,
                                             const std::vector<size_t> & variableSizes,
                                             const Parameters_ & params,
                                             const atlas::FieldSet & xb,
                                             const atlas::FieldSet & fg,
                                             const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistIncrOpSaberBlock starting" << std::endl;

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
  std::vector<std::string> requiredStateVariables{
    "potential_temperature",  // from file
    "exner",  // on theta levels from file ("exner_levels_minus_one" is on rho levels)
    "air_pressure",  // on theta levels from file ("air_pressure_levels_minus_one" is on rho levels)
    "air_temperature",  // to be populated in evalAirTemperature
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  // to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  // to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  // to be populated in evalSpecificHumidity
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in evalMassCloudLiquid
    "mass_content_of_cloud_ice_in_atmosphere_layer",  // to be populated in evalMassCloudIce
    "qrain",  // to be populated in evalMassRain
    "rht",  // to be populated in evalTotalRelativeHumidity
    "liquid_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "ice_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "cleff", "cfeff"  // to be populated in getMIOFields
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has_field(s)) {
      oops::Log::info() << "MoistIncrOpSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalMassCloudLiquid(augmentedStateFieldSet_);
  mo::evalMassCloudIce(augmentedStateFieldSet_);
  mo::evalMassRain(augmentedStateFieldSet_);
  mo::evalTotalRelativeHumidity(augmentedStateFieldSet_);
  mo::functions::getMIOFields(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::MoistIncrOpSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

MoistIncrOpSaberBlock::~MoistIncrOpSaberBlock() {
  oops::Log::trace() << classname() << "::~MoistIncrOpSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~MoistIncrOpSaberBlock");
  oops::Log::trace() << classname() << "::~MoistIncrOpSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("MoistIncrOpSaberBlock::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::qtTemperature2qqclqcfTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  mo::qqclqcf2qtTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::qtTemperature2qqclqcfAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  mo::qqclqcf2qtAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
