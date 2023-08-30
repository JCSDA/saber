/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydroBal.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/eval_air_pressure_levels.h"
#include "mo/eval_air_temperature.h"
#include "mo/eval_hydrostatic_balance.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/eval_virtual_potential_temperature.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "mo/model2geovals_varchange.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydroBal> makerHydroBal_("mo_hydro_bal");

// -----------------------------------------------------------------------------

HydroBal::HydroBal(const oops::GeometryData & outerGeometryData,
                   const oops::Variables & outerVars,
                   const eckit::Configuration & covarConf,
                   const Parameters_ & params,
                   const oops::FieldSet3D & xb,
                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydroBal starting" << std::endl;

  std::vector<std::string> requiredStateVariables{
    "air_temperature",
    "air_pressure",
    "air_pressure_levels_minus_one",
    "air_pressure_levels",
    "potential_temperature",   // from file
    "exner",  // from file on theta levels ("exner_levels_minus_one" is on rho levels)
    "exner_levels_minus_one",
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  //  to be populated in eval_sat_vapour_pressure_nl
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",
      //  to be populated in eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "virtual_potential_temperature",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels"
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.fieldSet().has(s)) {
      oops::Log::info() << "HydroBal variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }


  std::vector<std::string> requiredGeometryVariables{"height_levels"};
  for (const auto & s : requiredGeometryVariables) {
    if (outerGeometryData.fieldSet().has(s)) {
      augmentedStateFieldSet_.add(outerGeometryData.fieldSet()[s]);
    } else {
      augmentedStateFieldSet_.add(xb.fieldSet()[s]);
    }
  }

  // check how virtual potential temperature is calculated.
  mo::eval_air_pressure_levels_nl(augmentedStateFieldSet_);
  mo::eval_air_temperature_nl(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(params.svp_file, augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
              augmentedStateFieldSet_);
  mo::eval_virtual_potential_temperature_nl(augmentedStateFieldSet_);
  mo::evalHydrostaticExnerLevels(augmentedStateFieldSet_);
  mo::evalHydrostaticPressureLevels(augmentedStateFieldSet_);

  augmentedStateFieldSet_.haloExchange();

  oops::Log::trace() << classname() << "::HydroBal done" << std::endl;
}

// -----------------------------------------------------------------------------

HydroBal::~HydroBal() {
  oops::Log::trace() << classname() << "::~HydroBal starting" << std::endl;
  util::Timer timer(classname(), "~HydroBal");
  oops::Log::trace() << classname() << "::~HydroBal done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::eval_hydrobal_virtual_potential_temperature_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::eval_hydrobal_virtual_potential_temperature_ad(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  mo::eval_hydrobal_hydrostatic_exner_levels_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
