/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/GpToHp.h"

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
#include "mo/eval_exner.h"
#include "mo/eval_geostrophic_to_hydrostatic_pressure.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/eval_total_mixing_ratio.h"
#include "mo/eval_virtual_potential_temperature.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

namespace {

oops::Variables removeOuterOnlyVar(const oops::Variables & vars) {
  oops::Variables innerVars(vars);
  innerVars -= innerVars["hydrostatic_pressure_levels"];
  return innerVars;
}


}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GpToHp>
  makerGpToHp_("mo_hydrostatic_pressure_from_geostrophic_pressure");

// -----------------------------------------------------------------------------

GpToHp::GpToHp(const oops::GeometryData & outerGeometryData,
               const oops::Variables & outerVars,
               const eckit::Configuration & covarConf,
               const Parameters_ & params,
               const oops::FieldSet3D & xb,
               const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(removeOuterOnlyVar(getUnionOfInnerActiveAndOuterVars(params, outerVars))),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::GpToHp starting" << std::endl;
  // Covariance FieldSet
  covFieldSet_ = createGpRegressionStats(outerGeometryData.functionSpace(),
                                         outerGeometryData.fieldSet(),
                                         innerVars_,
                                         params.gptohpcovarianceparams.value());
  std::vector<std::string> requiredStateVariables{
    "air_temperature",
    "air_pressure_levels_minus_one",
    "exner_levels_minus_one",
    "exner",
    "potential_temperature",
    "air_pressure_levels",
    "air_pressure",
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in eval_total_mixing_ratio_nl
    "svp",  //  to be populated in eval_sat_vapour_pressure_nl
    "dlsvpdT",  //  to be populated in eval_derivative_ln_svp_wrt_temperature_nl
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",
      //  to be populated in eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "virtual_potential_temperature",
    "hydrostatic_exner_levels", "hydrostatic_pressure_levels"
     };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.fieldSet().has(s)) {
      oops::Log::info() << "HydrostaticExner variable " << s <<
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

  // we will need geometry here for height variables.
  mo::eval_air_pressure_levels_nl(augmentedStateFieldSet_);
  mo::eval_air_temperature_nl(augmentedStateFieldSet_);
  mo::eval_total_mixing_ratio_nl(augmentedStateFieldSet_);
  mo::eval_sat_vapour_pressure_nl(augmentedStateFieldSet_);
  mo::eval_derivative_ln_svp_wrt_temperature_nl(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
              augmentedStateFieldSet_);
  mo::eval_virtual_potential_temperature_nl(augmentedStateFieldSet_);
  mo::evalHydrostaticExnerLevels(augmentedStateFieldSet_);
  mo::evalHydrostaticPressureLevels(augmentedStateFieldSet_);
  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredCovarianceVariables;
  if (covFieldSet_.has("interpolation_weights")) {
    requiredCovarianceVariables.push_back("vertical_regression_matrices");
    requiredCovarianceVariables.push_back("interpolation_weights");
  }

  for (const auto & s : requiredCovarianceVariables) {
    augmentedStateFieldSet_.add(covFieldSet_[s]);
  }

  oops::Log::trace() << classname() << "::GpToHp done" << std::endl;
}

// -----------------------------------------------------------------------------

GpToHp::~GpToHp() {
  oops::Log::trace() << classname() << "::~GpToHp starting" << std::endl;
  util::Timer timer(classname(), "~GpToHp");
  oops::Log::trace() << classname() << "::~GpToHp done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_hydrostatic_pressure_levels_tl(fset.fieldSet(), augmentedStateFieldSet_);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_hydrostatic_pressure_levels_ad(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done"  << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHp::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  if (!fset.has("geostrophic_pressure_levels_minus_one")) {
    oops::Log::error() << "The inverse of block "
          << classname()
          << " is not correctly defined if geostrophic_pressure_levels_minus_one"
          << " is not provided as an input." << std::endl;
    throw eckit::UserError("Please only use leftInverseMultiply of this block "
                           "within the mo_hydrostatic_pressure block.", Here());
  }
  //   Allocate inner-only variables except air temperature
  oops::Variables innerOnlyVarsForInversion(innerOnlyVars_);
  innerOnlyVarsForInversion -= innerOnlyVarsForInversion["geostrophic_pressure_levels_minus_one"];
  checkFieldsAreNotAllocated(fset, innerOnlyVarsForInversion);
  allocateMissingFields(fset, innerOnlyVarsForInversion, innerOnlyVarsForInversion,
                        innerGeometryData_.functionSpace());

  // Retrieve unbalanced pressure from hydrostatic pressure and geostrophic pressure.
  mo::eval_hydrostatic_pressure_levels_tl_inv(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
