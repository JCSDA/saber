/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydrostaticExnerSaberBlock.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/HydrostaticExnerParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydrostaticExnerSaberBlock>
       makerHydrostaticExnerSaberBlock_("mo_hydrostatic_exner");

// -----------------------------------------------------------------------------

HydrostaticExnerSaberBlock::HydrostaticExnerSaberBlock(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydrostaticExnerSaberBlock starting" << std::endl;

  // Deserialize configuration
  HydrostaticExnerSaberBlockParameters params;
  params.deserialize(conf);

  // Input geometry and variables
  inputFunctionSpace_ = outputFunctionSpace;
  inputExtraFields_ = outputExtraFields;
  inputVars_ = params.outputVars.value();

  // Covariance FieldSet
  covFieldSet_ = createGpRegressionStats(outputFunctionSpace, outputExtraFields,
                                         activeVariableSizes, *params.activeVars.value(),
                                         params.hydrostaticexnerParams.value());

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

  std::vector<std::string> requiredGeometryVariables{"height_levels"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has_field(s)) {
      oops::Log::info() << "HydrostaticExnerSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(outputExtraFields[s]);
  }

  // we will need geometry here for height variables.
  mo::evalAirPressureLevels(augmentedStateFieldSet_);
  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);
  mo::evalHydrostaticExnerLevels(augmentedStateFieldSet_);
  mo::evalHydrostaticPressureLevels(augmentedStateFieldSet_);


  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredCovarianceVariables{
    "vertical_regression_matrices", "interpolation_weights"};

  for (const auto & s : requiredCovarianceVariables) {
    augmentedStateFieldSet_.add(covFieldSet_[s]);
  }

  oops::Log::trace() << classname() << "::HydrostaticExnerSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

HydrostaticExnerSaberBlock::~HydrostaticExnerSaberBlock() {
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~HydrostaticExnerSaberBlock");
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalHydrostaticPressureTL(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticExnerTL(fset, augmentedStateFieldSet_);
  const auto hydrostaticPressureView =
      atlas::array::make_view<const double, 2>(fset["hydrostatic_pressure_levels"]);
  auto airPressureView =
      atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  airPressureView.assign(hydrostaticPressureView);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  auto airPressureView =
      atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  auto hydrostaticPressureView =
      atlas::array::make_view<double, 2>(fset["hydrostatic_pressure_levels"]);

  for (atlas::idx_t jn = 0; jn < fset["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < fset["hydrostatic_exner_levels"].shape(1); ++jl) {
      hydrostaticPressureView(jn, jl) += airPressureView(jn, jl);
      airPressureView(jn, jl) = 0.0;
    }
  }
  mo::evalHydrostaticExnerAD(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticPressureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
