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

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
static SaberOuterBlockMaker<MoistIncrOpSaberBlock>
  makerMoistIncrOpSaberBlock_("mo_moistincrop");

// -----------------------------------------------------------------------------

MoistIncrOpSaberBlock::MoistIncrOpSaberBlock(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistIncrOpSaberBlock starting" << std::endl;

  // Deserialize configuration
  MoistIncrOpSaberBlockParameters params;
  params.deserialize(conf);

  // Input geometry and variables
  inputFunctionSpace_ = outputFunctionSpace;
  inputExtraFields_ = outputExtraFields;
  inputVars_ = params.outputVars.value();

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

void MoistIncrOpSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::qtTemperature2qqclqcfTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::qtTemperature2qqclqcfAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOpSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
