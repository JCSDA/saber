/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/MoistureControl.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/eval_air_temperature.h"
#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_cloud_liquid_mixing_ratio.h"
#include "mo/eval_moisture_control.h"
#include "mo/eval_moisture_incrementing_operator.h"
#include "mo/eval_rain_mixing_ratio.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/eval_total_relative_humidity.h"
#include "mo/eval_water_vapor_mixing_ratio.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/movader_covstats_interface.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<MoistureControl> makerMoistureControlBlock_("mo_moisture_control");

// -----------------------------------------------------------------------------

MoistureControl::MoistureControl(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::MoistureControl starting" << std::endl;

  // Covariance FieldSet
  covFieldSet_ = createMuStats(outerGeometryData.fieldSet(),
                               params.moistureControlParams.value());

  std::vector<std::string> requiredStateVariables{
    "air_temperature",
    "air_pressure",
    "potential_temperature",   // from file
    "exner",  // from file on theta levels ("exner_levels_minus_one" is on rho levels)
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in evalTotalMassMoistAir
    "svp",  //  to be populated in eval_sat_vapour_pressure_nl
    "dlsvpdT",  //  to be populated in eval_derivative_ln_svp_wrt_temperature_nl
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",
      //  to be populated in eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in
      // eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "mass_content_of_cloud_ice_in_atmosphere_layer",
      // to be populated in eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "qrain",  // to be populated in eval_rain_mixing_ratio_wrt_moist_air_and_condensed_water_nl
    "qt",  // to be populated in eval_total_water
    "rht",  // to be populated in eval_total_relative_humidity_nl
    "muA", "muH1",  // to be populated in function call from CovarianceStatisticsUtils.h
    "muRow1Column1", "muRow1Column2",  // to be populated in eval_moisture_control_traj
    "muRow2Column1", "muRow2Column2",  //   ""
    "muRecipDeterminant"  //   ""
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.fieldSet().has(s)) {
      oops::Log::info() << "MoistureControl variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  mo::eval_air_temperature_nl(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
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
  mo::eval_total_water(augmentedStateFieldSet_);
  mo::eval_total_relative_humidity_nl(augmentedStateFieldSet_);

  // populate "muA" and "muH1"
  interpMuStats(augmentedStateFieldSet_, covFieldSet_["muH1Stats"]);
  populateMuA(augmentedStateFieldSet_, covFieldSet_["muAStats"]);

  // populate "specific moisture control dependencies"
  mo::eval_moisture_control_traj(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::MoistureControl done" << std::endl;
}

// -----------------------------------------------------------------------------

MoistureControl::~MoistureControl() {
  oops::Log::trace() << classname() << "::~MoistureControl starting" << std::endl;
  util::Timer timer(classname(), "~MoistureControl");
  oops::Log::trace() << classname() << "::~MoistureControl done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_moisture_control_inv_tl(fset.fieldSet(), augmentedStateFieldSet_);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_moisture_control_inv_ad(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_moisture_control_tl(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const atlas::FieldSet & fields,
                              const MoistureControlCovarianceParameters & params) {
  // Get necessary parameters
  // path to covariance file with gp covariance parameters.
  std::string covFileName(params.covariance_file_path);
  // number of model levels
  std::size_t modelLevels;
  if (fields.has("height")) {
    modelLevels = fields["height"].shape(1);
  } else {
    modelLevels = fields["vert_coord"].shape(1);
  }

  // geostrophic pressure vertical regression statistics are grouped
  // into overlapping bins based on latitude;
  // number of bins associated with the gP vertical regression
  std::size_t muBins(static_cast<std::size_t>(params.mu_bins));

  // Need to setup derived state fields that we need.
  std::vector<std::string> shortnamesInFieldSet{
    "muAStats", "muH1Stats"};
  std::vector<std::string> shortnamesInFile{
    "M_inc_StdDev_binned", "H1_binned"};

  atlas::FieldSet statsFldSet;

  int sizeVec = static_cast<int>(modelLevels * muBins);
  std::vector<float> muStats1D(modelLevels * muBins, 0.0);

  // allocate and populate "muAStats", "muH1Stats"
  for (std::size_t i = 0; i < shortnamesInFile.size(); ++i) {
    covMuStats_f90(covFileName.size(),
                   covFileName.c_str(),
                   shortnamesInFile[i].size(),
                   shortnamesInFile[i].c_str(),
                   static_cast<int>(modelLevels),
                   muBins,
                   sizeVec,
                   muStats1D[0]);

    auto statsFld = atlas::Field(shortnamesInFieldSet[i],
      atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels, muBins));

    auto statsFldView = atlas::array::make_view<double, 2>(statsFld);
    std::size_t jn(0);
    for (std::size_t j = 0; j < modelLevels; ++j) {
      for (std::size_t b = 0; b < muBins; ++b, ++jn) {
        statsFldView(j, b) = static_cast<double>(muStats1D.at(jn));
      }
    }

    statsFldSet.add(statsFld);
  }

  return statsFldSet;
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
