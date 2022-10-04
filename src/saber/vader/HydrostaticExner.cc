/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydrostaticExner.h"

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
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydrostaticExner> makerHydrostaticExner_("mo_hydrostatic_exner");

// -----------------------------------------------------------------------------

HydrostaticExner::HydrostaticExner(const oops::GeometryData & outputGeometryData,
                                   const std::vector<size_t> & activeVariableSizes,
                                   const oops::Variables & outputVars,
                                   const Parameters_ & params,
                                   const atlas::FieldSet & xb,
                                   const atlas::FieldSet & fg,
                                   const std::vector<atlas::FieldSet> & fsetVec)
  : inputGeometryData_(outputGeometryData), inputVars_(outputVars), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydrostaticExner starting" << std::endl;

  // Get active variables
  oops::Variables activeVars = params.activeVars.value().get_value_or(outputVars);

  // Covariance FieldSet
  covFieldSet_ = createGpRegressionStats(outputGeometryData.functionSpace(),
                                         outputGeometryData.fieldSet(),
                                         activeVariableSizes,
                                         activeVars,
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
      oops::Log::info() << "HydrostaticExner variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(outputGeometryData.fieldSet()[s]);
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

  oops::Log::trace() << classname() << "::HydrostaticExner done" << std::endl;
}

// -----------------------------------------------------------------------------

HydrostaticExner::~HydrostaticExner() {
  oops::Log::trace() << classname() << "::~HydrostaticExner starting" << std::endl;
  util::Timer timer(classname(), "~HydrostaticExner");
  oops::Log::trace() << classname() << "::~HydrostaticExner done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExner::multiply(atlas::FieldSet & fset) const {
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

void HydrostaticExner::multiplyAD(atlas::FieldSet & fset) const {
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

void HydrostaticExner::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExner::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

/// \details This extracts the hp_gp regression matrix for a number of
///          of overlapping latitude bands from the operational covariance
///          statistics file. The matrices are stored as a single field.
///
/// B = (vertical regression matrix bin_0)
///     (vertical regression matrix bin_1)
///     (          ...                   )
///     (vertical regression matrix bin_m)
/// Since each matrix is square we can easily infer the bin index from the row index
/// First index of vertRegView is bin_index * number of levels + level index,
///     the second is number of levels associated with matrix column.
///
/// The interpolation weights are calculated for each grid point location
/// The second index relates to the bin index.
/// We ensure that across all bins for a grid point we sum to 1.
///
atlas::FieldSet createGpRegressionStats(const atlas::FunctionSpace & functionSpace,
                                        const atlas::FieldSet & extraFields,
                                        const std::vector<size_t> & variableSizes,
                                        const oops::Variables & inputVars,
                                        const HydrostaticExnerCovarianceParameters & params) {
  // Get necessary parameters
  // path to covariance file with gp covariance parameters.
  std::string covFileName(params.covariance_file_path);
  // number of latitudes that existed in the generation of the covariance file
  std::size_t covGlobalNLats(static_cast<std::size_t>(params.covariance_nlat));
  // number of model levels
  std::size_t modelLevels(variableSizes[0]);
  // geostrophic pressure vertical regression statistics are grouped
  // into overlapping bins based on latitude;
  // number of bins associated with the gP vertical regression
  std::size_t gPBins(static_cast<std::size_t>(params.gp_regression_bins));

  atlas::FieldSet gpStatistics;

  gpStatistics.add(createGpRegressionMatrices(covFileName, gPBins, modelLevels));

  gpStatistics.add(createGpRegressionWeights(functionSpace, extraFields,
                                             covFileName, covGlobalNLats, gPBins));

  return gpStatistics;
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
