/*
 * (C) Crown Copyright 2022-2024 Met Office
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

#include "mo/common_varchange.h"
#include "mo/eval_air_temperature.h"
#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_cloud_liquid_mixing_ratio.h"
#include "mo/eval_mio_fields.h"
#include "mo/eval_moisture_incrementing_operator.h"
#include "mo/eval_rain_mixing_ratio.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/eval_total_mixing_ratio.h"
#include "mo/eval_total_relative_humidity.h"
#include "mo/eval_water_vapor_mixing_ratio.h"
#include "mo/functions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<MoistIncrOp> makerMoistIncrOp_("mo_moistincrop");

// -----------------------------------------------------------------------------

MoistIncrOp::MoistIncrOp(const oops::GeometryData & outerGeometryData,
                         const oops::Variables & outerVars,
                         const eckit::Configuration & covarConf,
                         const Parameters_ & params,
                         const oops::FieldSet3D & xb,
                         const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistIncrOp starting" << std::endl;

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredStateVariables{
    "potential_temperature",  // from file
    "exner",  // on theta levels from file ("exner_levels_minus_one" is on rho levels)
    "air_pressure",  // on theta levels from file ("air_pressure_levels_minus_one" is on rho levels)
    "air_temperature",  // to be populated in eval_air_temperature_nl
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  // to be populated in eval_total_mixing_ratio_nl
    "svp",  //  to be populated in eval_sat_vapour_pressure_nl
    "dlsvpdT",  //  to be populated in eval_derivative_ln_svp_wrt_temperature_nl
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",
      // to be populated in eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in
      // eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "mass_content_of_cloud_ice_in_atmosphere_layer",
      // to be populated in eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "qrain",  // to be populated in eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "rht",  // to be populated in eval_total_relative_humidity_nl
    "liquid_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "ice_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "cleff", "cfeff"  // to be populated in eval_mio_fields_nl
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.fieldSet().has(s)) {
      oops::Log::info() << "MoistIncrOp variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  mo::eval_air_temperature_nl(augmentedStateFieldSet_);
  mo::eval_total_mixing_ratio_nl(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(augmentedStateFieldSet_);
  mo::eval_derivative_ln_svp_wrt_temperature_nl(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
              augmentedStateFieldSet_);
  mo::eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
              augmentedStateFieldSet_);
  mo::eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
              augmentedStateFieldSet_);
  mo::eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_nl(augmentedStateFieldSet_);
  mo::eval_total_relative_humidity_nl(augmentedStateFieldSet_);
  mo::eval_mio_fields_nl(params.mio_file, augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::MoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

MoistIncrOp::~MoistIncrOp() {
  oops::Log::trace() << classname() << "::~MoistIncrOp starting" << std::endl;
  util::Timer timer(classname(), "~MoistIncrOp");
  oops::Log::trace() << classname() << "::~MoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_moisture_incrementing_operator_tl(fset.fieldSet(), augmentedStateFieldSet_);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_moisture_incrementing_operator_ad(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  if (!fset.has("air_temperature")) {
    oops::Log::error() << "The inverse of the moisture incrementing operator "
          << "is not correctly defined if air_temperature is not provided "
          << "as an input." << std::endl;
    throw eckit::UserError("Please only use leftInverseMultiply of the mo_moistincrop block "
                           "within the mo_super_mio block.", Here());
  }
  //   Allocate inner-only variables except air temperature
  oops::Variables innerOnlyVarsForInversion(innerOnlyVars_);
  innerOnlyVarsForInversion -= "air_temperature";
  checkFieldsAreNotAllocated(fset, innerOnlyVarsForInversion);
  allocateMissingFields(fset, innerOnlyVarsForInversion, innerOnlyVarsForInversion,
                        innerGeometryData_.functionSpace());

  mo::eval_total_water_tl(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistIncrOp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
