/*
 * (C) Crown Copyright 2022 Met Office
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
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/eval_sat_vapour_pressure.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/movader_covstats_interface.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<MoistureControl> makerMoistureControlBlock_("mo_moisture_control");

// -----------------------------------------------------------------------------

MoistureControl::MoistureControl(const oops::GeometryData & outerGeometryData,
                                 const std::vector<size_t> & activeVariableSizes,
                                 const oops::Variables & outerVars,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const std::vector<atlas::FieldSet> & fsetVec)
  : innerGeometryData_(outerGeometryData), innerVars_(outerVars), augmentedStateFieldSet_()
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
    "svp", "dlsvpdT",  //  to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  //  to be populated in evalSpecificHumidity
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in evalMassCloudLiquid
    "mass_content_of_cloud_ice_in_atmosphere_layer",  // to be populated in evalMassCloudIce
    "qrain",  // to be populated in evalMassRain
    "qt",  // to be populated in qqclqcf2qtTL
    "rht",  // to be populated in evalTotalRelativeHumidity
    "muA", "muH1",  // to be populated in function call from CovarianceStatisticsUtils.h
    "muRow1Column1", "muRow1Column2",  // to be populated in evalMoistureControlDependencies
    "muRow2Column1", "muRow2Column2",  //   ""
    "muRecipDeterminant"  //   ""
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.has(s)) {
      oops::Log::info() << "MoistureControl variable " << s <<
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
  mo::qqclqcf2qt(augmentedStateFieldSet_);
  mo::evalTotalRelativeHumidity(augmentedStateFieldSet_);

  // populate "muA" and "muH1"
  for (auto & covFld : covFieldSet_) {
    populateInterpMuStats(augmentedStateFieldSet_,
                          covFld);
  }

  // populate "specific moisture control dependencies"
  mo::evalMoistureControlDependencies(augmentedStateFieldSet_);

  augmentedStateFieldSet_.haloExchange();

  oops::Log::trace() << classname() << "::MoistureControl done" << std::endl;
}

// -----------------------------------------------------------------------------

MoistureControl::~MoistureControl() {
  oops::Log::trace() << classname() << "::~MoistureControl starting" << std::endl;
  util::Timer timer(classname(), "~MoistureControl");
  oops::Log::trace() << classname() << "::~MoistureControl done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalQtThetaTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalQtThetaAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const atlas::FieldSet & extraFields,
                              const MoistureControlCovarianceParameters & params) {
  // Get necessary parameters
  // path to covariance file with gp covariance parameters.
  std::string covFileName(params.covariance_file_path);
  // number of model levels
  std::size_t modelLevels;
  if (extraFields.has("height")) {
    modelLevels = extraFields["height"].levels();
  } else {
    modelLevels = extraFields["vunit"].levels();
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
