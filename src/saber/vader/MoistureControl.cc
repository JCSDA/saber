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

#include "mo/constants.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/movader_covstats_interface.h"

namespace saber {
namespace vader {

namespace {

using atlas::array::make_view;
using atlas::idx_t;

void eval_moisture_control_inv_tl(atlas::FieldSet & incFlds,
                                  const atlas::FieldSet & augStateFlds) {
  // Using Cramer's rule to calculate inverse.
  auto muRecipDeterView = make_view<const double, 2>(augStateFlds["muRecipDeterminant"]);
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View  = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto muIncView = make_view<const double, 2>(incFlds["mu"]);
  auto thetavIncView = make_view<const double, 2>(incFlds["virtual_potential_temperature"]);
  auto qtIncView = make_view<double, 2>(incFlds["qt"]);
  auto thetaIncView = make_view<double, 2>(incFlds["air_potential_temperature"]);

  const atlas::idx_t n_levels(incFlds["mu"].shape(1));
  const idx_t sizeOwned =
        util::getSizeOwned(incFlds["mu"].functionspace());
  atlas_omp_parallel_for(idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      // VAR equivalent in Var_UpPFtheta_qT.f90 for thetaIncView
      // (beta2 * muA * theta_v' +   beta1 * mu') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      thetaIncView(ih, ilev) = muRecipDeterView(ih, ilev) * (
                             muRow1Column1View(ih, ilev) * thetavIncView(ih, ilev)
                           - muRow2Column1View(ih, ilev) * muIncView(ih, ilev) );

      // VAR equivalent in Var_UpPFtheta_qT.f90 for qtIncView
      // (alpha1 * mu_v' -   alpha2 * muA * thetav') /
      // (alpha1 * beta2 * muA - alpha2 * muA * beta1)
      qtIncView(ih, ilev) = muRecipDeterView(ih, ilev) * (
                           muRow2Column2View(ih, ilev) * muIncView(ih, ilev) -
                           muRow1Column2View(ih, ilev) * thetavIncView(ih, ilev) );
    }
  }
  incFlds["air_potential_temperature"].set_dirty();
  incFlds["qt"].set_dirty();
}

void eval_moisture_control_inv_ad(atlas::FieldSet & hatFlds,
                                  const atlas::FieldSet & augStateFlds) {
  auto muRecipDeterView = make_view<const double, 2>(augStateFlds["muRecipDeterminant"]);
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View  = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto qtHatView = make_view<double, 2>(hatFlds["qt"]);
  auto muHatView = make_view<double, 2>(hatFlds["mu"]);
  auto thetavHatView = make_view<double, 2>(hatFlds["virtual_potential_temperature"]);
  auto thetaHatView = make_view<double, 2>(hatFlds["air_potential_temperature"]);

  const atlas::idx_t n_levels(hatFlds["mu"].shape(1));
  const idx_t sizeOwned =
        util::getSizeOwned(hatFlds["mu"].functionspace());
  atlas_omp_parallel_for(idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      thetavHatView(ih, ilev) += muRecipDeterView(ih, ilev) *
                                 muRow1Column1View(ih, ilev) * thetaHatView(ih, ilev);
      muHatView(ih, ilev) -= muRecipDeterView(ih, ilev) *
                             muRow2Column1View(ih, ilev) * thetaHatView(ih, ilev);
      thetavHatView(ih, ilev) -= muRecipDeterView(ih, ilev) *
                                 muRow1Column2View(ih, ilev) * qtHatView(ih, ilev);
      muHatView(ih, ilev) += muRecipDeterView(ih, ilev) *
                             muRow2Column2View(ih, ilev) * qtHatView(ih, ilev);
      thetaHatView(ih, ilev) = 0.0;
      qtHatView(ih, ilev) = 0.0;
    }
  }
  hatFlds["air_potential_temperature"].set_dirty();
  hatFlds["qt"].set_dirty();
  hatFlds["mu"].set_dirty();
  hatFlds["virtual_potential_temperature"].set_dirty();
}


void eval_moisture_control_tl(atlas::FieldSet & incFlds,
                              const atlas::FieldSet & augStateFlds) {
  auto muRow1Column1View = make_view<const double, 2>(augStateFlds["muRow1Column1"]);
  auto muRow1Column2View = make_view<const double, 2>(augStateFlds["muRow1Column2"]);
  auto muRow2Column1View = make_view<const double, 2>(augStateFlds["muRow2Column1"]);
  auto muRow2Column2View = make_view<const double, 2>(augStateFlds["muRow2Column2"]);
  auto thetaIncView = make_view<const double, 2>(incFlds["air_potential_temperature"]);
  auto qtIncView = make_view<const double, 2>(incFlds["qt"]);
  auto muIncView = make_view<double, 2>(incFlds["mu"]);
  auto thetavIncView = make_view<double, 2>(incFlds["virtual_potential_temperature"]);

  const idx_t n_levels(incFlds["mu"].shape(1));
  const idx_t sizeOwned =
        util::getSizeOwned(incFlds["mu"].functionspace());
  atlas_omp_parallel_for(idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      muIncView(ih, ilev) = muRow1Column1View(ih, ilev) * qtIncView(ih, ilev)
                          + muRow1Column2View(ih, ilev) * thetaIncView(ih, ilev);
      thetavIncView(ih, ilev) = muRow2Column1View(ih, ilev) * qtIncView(ih, ilev)
                            + muRow2Column2View(ih, ilev) * thetaIncView(ih, ilev);
    }
  }
  incFlds["mu"].set_dirty();
  incFlds["virtual_potential_temperature"].set_dirty();
}

void eval_moisture_control_traj(atlas::FieldSet & fields) {
  auto qtView = make_view<const double, 2>(fields["qt"]);
  auto qView = make_view<const double, 2>(
                              fields["water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"]);
  auto thetaView = make_view<const double, 2>(fields["air_potential_temperature"]);
  auto exnerView = make_view<const double, 2>(fields["dimensionless_exner_function"]);
  auto dlsvpdTView = make_view<const double, 2>(fields["dlsvpdT"]);
  auto qsatView = make_view<const double, 2>(fields["qsat"]);
  auto muAView = make_view<const double, 2>(fields["muA"]);
  auto muH1View = make_view<const double, 2>(fields["muH1"]);

  // muRow1Column1, muRow1Column2, muRow2Column1, muRow2Column2
  // are coefficients of a (2x2) matrix = A
  //  (mu')       = A (qt')     where A is
  //  (theta_v')      (theta')
  //
  //  ( muA/qsat    - (muA/qsat) muH1 qT exner_bar dlsvpdT )
  //  (                                                    )
  //  (c_v theta     (1 + cv q)                            )
  //
  auto muRow1Column1View = make_view<double, 2>(fields["muRow1Column1"]);
  auto muRow1Column2View = make_view<double, 2>(fields["muRow1Column2"]);
  auto muRow2Column1View = make_view<double, 2>(fields["muRow2Column1"]);
  auto muRow2Column2View = make_view<double, 2>(fields["muRow2Column2"]);
  auto muRecipDeterminantView = make_view<double, 2>(fields["muRecipDeterminant"]);

  // the comments below are there to allow checking with the VAR code.
  const idx_t n_levels(fields["air_potential_temperature"].shape(1));
  const idx_t sizeOwned =
        util::getSizeOwned(fields["air_potential_temperature"].functionspace());
  atlas_omp_parallel_for(idx_t ih = 0; ih < sizeOwned; ++ih) {
    for (idx_t ilev = 0; ilev < n_levels; ++ilev) {
      muRow1Column1View(ih, ilev) = muAView(ih, ilev) / qsatView(ih, ilev);  // beta2 * muA
      muRow1Column2View(ih, ilev) = - qtView(ih, ilev)  * muH1View(ih, ilev)
        * exnerView(ih, ilev) * dlsvpdTView(ih, ilev) * muRow1Column1View(ih, ilev);
      // alpha2 * muA
      muRow2Column1View(ih, ilev) = ::mo::constants::c_virtual * thetaView(ih, ilev);   // beta1
      muRow2Column2View(ih, ilev) = 1.0 + ::mo::constants::c_virtual * qView(ih, ilev);  // alpha1
      muRecipDeterminantView(ih, ilev) = 1.0 /(
        muRow2Column2View(ih, ilev) * muRow1Column1View(ih, ilev)
        - muRow1Column2View(ih, ilev) * muRow2Column1View(ih, ilev));
           // 1/( alpha1 * beta2 * muA - alpha2 * muA * beta1)
    }
  }
  fields["muRow1Column1"].set_dirty();
  fields["muRow1Column2"].set_dirty();
  fields["muRow2Column1"].set_dirty();
  fields["muRow2Column2"].set_dirty();
  fields["muRecipDeterminant"].set_dirty();
}

}  // namespace

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
    nlevs_(xb["air_temperature"].levels()),
    params_(params),
    covFieldSet_(),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistureControl starting" << std::endl;

  // copy all required variables from the background fieldset
  const oops::Variables mandatoryStateVariables = params.mandatoryStateVars();
  augmentedStateFieldSet_.clear();
  for (const auto & s : mandatoryStateVariables.variables()) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

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
  eval_moisture_control_inv_tl(fset.fieldSet(), augmentedStateFieldSet_);

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

  eval_moisture_control_inv_ad(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_moisture_control_tl(fset.fieldSet(), augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::read() {
  oops::Log::trace() << classname() << "::read start " <<  std::endl;
  MoistureControlReadParameters mparams;
  const auto & calibparams = params_.calibrationParams.value();
  if (calibparams != boost::none) {
    const auto & calibrationReadParams = calibparams->calibrationReadParams.value();
    if (calibrationReadParams != boost::none) {
       mparams = calibrationReadParams.value();
    }
  } else {
    mparams = *params_.readParams.value();
  }

  // Covariance FieldSet
  covFieldSet_ = createMuStats(nlevs_,
                               innerGeometryData_.fieldSet(),
                               mparams);

  std::vector<std::string> additionalStateVariables{
    "muA", "muH1",  // to be populated in function call from CovarianceStatisticsUtils.h
    "muRow1Column1", "muRow1Column2",  // to be populated in eval_moisture_control_traj
    "muRow2Column1", "muRow2Column2",  //   ""
    "muRecipDeterminant"  //   ""
  };

  // create fields for temporary variables required here
  for (const auto & s : additionalStateVariables) {
    atlas::Field field = innerGeometryData_.functionSpace()->createField<double>(
        atlas::option::name(s) | atlas::option::levels(nlevs_));
    augmentedStateFieldSet_.add(field);
  }

  // populate "muA" and "muH1"
  interpMuStats(augmentedStateFieldSet_, covFieldSet_["muH1Stats"]);
  populateMuA(augmentedStateFieldSet_, covFieldSet_["muAStats"]);
  // populate "specific moisture control dependencies"
  eval_moisture_control_traj(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::directCalibration(const oops::FieldSets & fset) {
  oops::Log::trace() << classname() << "::directCalibration start" << std::endl;
  const auto & calibparams = params_.calibrationParams.value();
  ASSERT(calibparams != boost::none);
  const auto & calibrationReadParams = calibparams->calibrationReadParams.value();
  if (calibrationReadParams != boost::none) {
    MoistureControl::read();
  }
  oops::Log::trace() << classname() << "::directCalibration end" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControl::write() const {
  oops::Log::trace() << classname() << "::write start" << std::endl;

  oops::Log::trace() << classname() << "::write end" << std::endl;
}


// -----------------------------------------------------------------------------

void MoistureControl::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const size_t & modelLevelsDefault,
                              const atlas::FieldSet & fields,
                              const MoistureControlReadParameters & params) {
  // Get necessary parameters
  // path to covariance file with gp covariance parameters.
  std::string covFileName(params.covariance_file_path);
  // number of model levels
  std::size_t modelLevels;
  if (fields.has("height_above_mean_sea_level")) {
    modelLevels = fields["height_above_mean_sea_level"].shape(1);
  } else {
    modelLevels = modelLevelsDefault;
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
