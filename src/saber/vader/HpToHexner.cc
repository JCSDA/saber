/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HpToHexner.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/eval_exner.h"
#include "mo/eval_geostrophic_to_hydrostatic_pressure.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HpToHexner>
  makerHpToHexner_("mo_hydrostatic_pressure_to_hydrostatic_exner");

// -----------------------------------------------------------------------------

HpToHexner::HpToHexner(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const atlas::FieldSet & xb,
                       const atlas::FieldSet & fg,
                       const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HpToHexner starting" << std::endl;

  std::vector<std::string> requiredStateVariables{
    "air_temperature",
    "air_pressure_levels_minus_one",
    "exner_levels_minus_one",
    "exner",
    "potential_temperature",
    "air_pressure_levels",
    "air_pressure",
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  //  to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  //  to be populated in evalSpecificHumidity
    "virtual_potential_temperature",
    "hydrostatic_exner_levels", "hydrostatic_pressure_levels"
     };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has(s)) {
      oops::Log::info() << "HpToHexner variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }


  std::vector<std::string> requiredGeometryVariables{"height_levels"};
  for (const auto & s : requiredGeometryVariables) {
    if (outerGeometryData.fieldSet().has(s)) {
      augmentedStateFieldSet_.add(outerGeometryData.fieldSet()[s]);
    } else {
      augmentedStateFieldSet_.add(xb[s]);
    }
  }

  // we will need geometry here for height variables.
  mo::evalAirPressureLevels(augmentedStateFieldSet_);
  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(params.svp_file, augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);
  mo::evalHydrostaticExnerLevels(augmentedStateFieldSet_);
  mo::evalHydrostaticPressureLevels(augmentedStateFieldSet_);

  augmentedStateFieldSet_.haloExchange();

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredCovarianceVariables;
  if (covFieldSet_.has("interpolation_weights")) {
    requiredCovarianceVariables.push_back("vertical_regression_matrices");
    requiredCovarianceVariables.push_back("interpolation_weights");
  }

  for (const auto & s : requiredCovarianceVariables) {
    augmentedStateFieldSet_.add(covFieldSet_[s]);
  }

  oops::Log::trace() << classname() << "::HpToHexner done" << std::endl;
}

// -----------------------------------------------------------------------------

HpToHexner::~HpToHexner() {
  oops::Log::trace() << classname() << "::~HpToHexner starting" << std::endl;
  util::Timer timer(classname(), "~HpToHexner");
  oops::Log::trace() << classname() << "::~HpToHexner done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::eval_hydrostatic_exner_levels_tl(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::eval_hydrostatic_exner_levels_ad(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Retrieve hydrostatic pressure from hydrostatic Exner.
  mo::eval_hydrostatic_exner_levels_tl_inv(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
