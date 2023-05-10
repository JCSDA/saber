/*
 * (C) Crown Copyright 2022-2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/MoistIncrOp.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "../../vader/src/mo/common_varchange.h"
#include "../../vader/src/mo/control2analysis_linearvarchange.h"
#include "../../vader/src/mo/eval_mio_fields.h"
#include "../../vader/src/mo/eval_moisture_incrementing_operator.h"
#include "../../vader/src/mo/eval_sat_vapour_pressure.h"
#include "../../vader/src/mo/eval_total_relative_humidity.h"
#include "../../vader/src/mo/functions.h"
#include "../../vader/src/mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<MoistIncrOp> makerMoistIncrOp_("mo_moistincrop");

// -----------------------------------------------------------------------------

MoistIncrOp::MoistIncrOp(const oops::GeometryData & outerGeometryData,
                         const std::vector<size_t> & activeVariableSizes,
                         const oops::Variables & outerVars,
                         const eckit::Configuration & covarConf,
                         const Parameters_ & params,
                         const atlas::FieldSet & xb,
                         const atlas::FieldSet & fg,
                         const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistIncrOp starting" << std::endl;

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
    "rht",  // to be populated in eval_total_relative_humidity_nl
    "liquid_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "ice_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "cleff", "cfeff"  // to be populated in getMIOFields
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has(s)) {
      oops::Log::info() << "MoistIncrOp variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(params.svp_file, augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalMassCloudLiquid(augmentedStateFieldSet_);
  mo::evalMassCloudIce(augmentedStateFieldSet_);
  mo::evalMassRain(augmentedStateFieldSet_);
  mo::eval_total_relative_humidity_nl(augmentedStateFieldSet_);
  mo::eval_mio_fields_nl(params.mio_file, augmentedStateFieldSet_);

  augmentedStateFieldSet_.haloExchange();

  oops::Log::trace() << classname() << "::MoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

MoistIncrOp::~MoistIncrOp() {
  oops::Log::trace() << classname() << "::~MoistIncrOp starting" << std::endl;
  util::Timer timer(classname(), "~MoistIncrOp");
  oops::Log::trace() << classname() << "::~MoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::eval_moisture_incrementing_operator_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::eval_moisture_incrementing_operator_ad(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  mo::eval_total_water_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
