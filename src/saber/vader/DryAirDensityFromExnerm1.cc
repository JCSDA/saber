/*
 * (C) Crown Copyright 2024-2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensityFromExnerm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "mo/constants.h"
#include "mo/eval_dry_air_density.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/for_each.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"


namespace saber {
namespace vader {

namespace {

using View = atlas::array::LocalView<double, 1>;
using ConstView = atlas::array::LocalView<const double, 1>;

const char specific_humidity_mo[] = "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water";

/// \details Calculate the dry air density from potential temperature,
///          specific humidity and exner pressure (using exner_pressure_levels_minus_one)
void eval_dry_air_density_from_exner_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                         const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_pressure_levels_minus_one_tl()] starting ..."
                     << std::endl;

  util::for_each_column(
      [=](ConstView dryrho,
          ConstView exner,
          ConstView hl,
          ConstView h,
          ConstView pt,
          ConstView q,
          ConstView qcl,
          ConstView qcf,
          ConstView exnerInc,
          ConstView ptInc,
          ConstView qInc,
          ConstView qclInc,
          ConstView qcfInc,
          View dryrhoInc) {
        double h_minus_hl;
        double hl_minus_hm1;
        double vptdrydens;
        double vptdrydensInc;
        double vptdrydens_jlm1;
        double vptdrydensInc_jlm1;
        double vptdrydens_intp_times_h_minus_hm1;
        double vptdrydensInc_intp_times_h_minus_hm1;

        vptdrydens = pt(0) * (1.0 + ::mo::constants::c_virtual * q(0)
                     - qcl(0) - qcf(0)) / (1.0 - q(0)
                     - qcl(0) - qcf(0));
        vptdrydensInc = ((1.0 + ::mo::constants::c_virtual * q(0)
                          - qcl(0) - qcf(0)) * ptInc(0)
                         + ((1.0 + ::mo::constants::c_virtual) * qInc(0)
                            * (1.0 - qcl(0) - qcf(0)) * pt(0)
                            / (1.0 - q(0) - qcl(0) - qcf(0)))
                         + (qclInc(0) *  q(0)
                           * (1.0 + ::mo::constants::c_virtual) * pt(0))
                           / (1.0 - q(0) - qcl(0) - qcf(0))
                         + (qcfInc(0) *  q(0)
                           * (1.0 + ::mo::constants::c_virtual) * pt(0))
                           / (1.0 - q(0) - qcl(0) - qcf(0)))
                        / (1.0 - q(0) - qcl(0) - qcf(0));

        dryrhoInc(0) = dryrho(0) * (
          (1.0 - ::mo::constants::rd_over_cp)  * exnerInc(0) /
          (exner(0) * ::mo::constants::rd_over_cp) -
          vptdrydensInc / vptdrydens);

        const atlas::idx_t nb_levels = dryrho.shape(0);  // in this scope, dryrho is a column view
        for (atlas::idx_t jl = 1; jl < nb_levels; ++jl) {
          h_minus_hl = h(jl) - hl(jl);
          hl_minus_hm1 = hl(jl) - h(jl-1);
          vptdrydens_jlm1 = vptdrydens;
          vptdrydens = pt(jl) * (1.0 + ::mo::constants::c_virtual * q(jl)
                       - qcl(jl) - qcf(jl)) / (1.0 - q(jl)
                       - qcl(jl) - qcf(jl));
          vptdrydensInc_jlm1 = vptdrydensInc;
          vptdrydensInc = ((1.0 + ::mo::constants::c_virtual * q(jl)
                            - qcl(jl) - qcf(jl)) * ptInc(jl)
                           + ((1.0 + ::mo::constants::c_virtual) * qInc(jl)
                              * (1.0 - qcl(jl) - qcf(jl)) * pt(jl)
                              / (1.0 - q(jl) - qcl(jl) - qcf(jl)))
                           + (qclInc(jl) *  q(jl)
                              * (1.0 + ::mo::constants::c_virtual) * pt(jl))
                              / (1.0 - q(jl) - qcl(jl) - qcf(jl))
                           + (qcfInc(jl) *  q(jl)
                              * (1.0 + ::mo::constants::c_virtual) * pt(jl))
                              / (1.0 - q(jl) - qcl(jl) - qcf(jl)))
                          / (1.0 - q(jl) - qcl(jl) - qcf(jl));
          vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                               hl_minus_hm1 * vptdrydens;
          vptdrydensInc_intp_times_h_minus_hm1 = hl_minus_hm1 * vptdrydensInc +
                                                 h_minus_hl * vptdrydensInc_jlm1;

          dryrhoInc(jl) = dryrho(jl) * (
            (1.0 - ::mo::constants::rd_over_cp) *
            (exnerInc(jl) /
            (exner(jl) * ::mo::constants::rd_over_cp)) -
            vptdrydensInc_intp_times_h_minus_hm1 / vptdrydens_intp_times_h_minus_hm1);
        }
      },
      stateFlds["dry_air_density_levels_minus_one"],
      stateFlds["dimensionless_exner_function_levels_minus_one"],
      stateFlds["height_above_mean_sea_level_levels"],
      stateFlds["height_above_mean_sea_level"],
      stateFlds["air_potential_temperature"],
      stateFlds[specific_humidity_mo],
      stateFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"],
      stateFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"],
      incFlds["dimensionless_exner_function_levels_minus_one"],
      incFlds["air_potential_temperature"],
      incFlds[specific_humidity_mo],
      incFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"],
      incFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"],
      incFlds["dry_air_density_levels_minus_one"]);

  incFlds["dry_air_density_levels_minus_one"].set_dirty();

  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_tl()] ... exit"
                     << std::endl;
}

// -------------------------------------------------------------------------------------------------

void eval_dry_air_density_from_exner_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                            const atlas::FieldSet & stateFlds) {
  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_ad()] starting ..."
                     << std::endl;


  util::for_each_column(
      [=](ConstView dryrho,
          ConstView exner,
          ConstView hl,
          ConstView h,
          ConstView pt,
          ConstView q,
          ConstView qcl,
          ConstView qcf,
          View exnerHat,
          View ptHat,
          View qHat,
          View qclHat,
          View qcfHat,
          View dryrhoHat) {
        double h_minus_hl;
        double hl_minus_hm1;
        double vptdrydens;
        double vptdrydens_intp_times_h_minus_hm1;
        double vptdrydensHat;
        double vptdrydens_jlm1;
        double vptdrydensHat_jlm1;

        const atlas::idx_t nb_levels = dryrho.shape(0);
        for (atlas::idx_t jl = nb_levels-1; jl >= 1; --jl) {
          // Passive fields.
          h_minus_hl = h(jl) - hl(jl);
          hl_minus_hm1 = hl(jl) - h(jl-1);
          vptdrydens = pt(jl) * (1.0 + ::mo::constants::c_virtual * q(jl)
                       - qcl(jl) - qcf(jl)) / (1.0 - q(jl)
                       - qcl(jl) - qcf(jl));

          exnerHat(jl) += dryrho(jl) *
            (1.0 - ::mo::constants::rd_over_cp)  * dryrhoHat(jl) /
            (exner(jl) * ::mo::constants::rd_over_cp);

          vptdrydens_jlm1 = pt(jl-1) * (1.0 + ::mo::constants::c_virtual * q(jl-1)
                            - qcl(jl-1) - qcf(jl-1)) / (1.0 - q(jl-1)
                            - qcl(jl-1) - qcf(jl-1));
          vptdrydens_intp_times_h_minus_hm1 = h_minus_hl * vptdrydens_jlm1 +
                                               hl_minus_hm1 * vptdrydens;
          vptdrydensHat = - dryrho(jl) * dryrhoHat(jl) * hl_minus_hm1
                          / vptdrydens_intp_times_h_minus_hm1;
          vptdrydensHat_jlm1 = - dryrho(jl) * dryrhoHat(jl) * h_minus_hl /
                               vptdrydens_intp_times_h_minus_hm1;
          ptHat(jl) += (1.0 + ::mo::constants::c_virtual * q(jl)
                                - qcl(jl) - qcf(jl))
                               / (1.0 - q(jl) - qcl(jl) - qcf(jl))
                               * vptdrydensHat;
          qHat(jl) += (1.0 + ::mo::constants::c_virtual) * pt(jl)
                               * (1.0 - qcl(jl) - qcf(jl))
                               / ((1.0 - q(jl) - qcl(jl) - qcf(jl))
                                   * (1.0 - q(jl) - qcl(jl) - qcf(jl)))
                               * vptdrydensHat;
          qclHat(jl) += pt(jl) *  q(jl) * (1.0 + ::mo::constants::c_virtual)
                                / ((1.0 - q(jl) - qcl(jl) - qcf(jl))
                                    * (1.0 - q(jl) - qcl(jl) - qcf(jl)))
                                * vptdrydensHat;
          qcfHat(jl) += pt(jl) *  q(jl) * (1.0 + ::mo::constants::c_virtual)
                                / ((1.0 - q(jl) - qcl(jl) - qcf(jl))
                                    * (1.0 - q(jl) - qcl(jl) - qcf(jl)))
                                * vptdrydensHat;
          ptHat(jl-1) += (1.0 + ::mo::constants::c_virtual * q(jl-1)
                                  - qcl(jl-1) - qcf(jl-1))
                                 / (1.0 - q(jl-1) - qcl(jl-1) - qcf(jl-1))
                                 * vptdrydensHat_jlm1;
          qHat(jl-1) += (1.0 + ::mo::constants::c_virtual) * pt(jl-1)
                                * (1.0 - qcl(jl-1) - qcf(jl-1))
                                / ((1.0 - q(jl-1) - qcl(jl-1) - qcf(jl-1))
                                    * (1.0 - q(jl-1) - qcl(jl-1) - qcf(jl-1)))
                                * vptdrydensHat_jlm1;
          qclHat(jl-1) += pt(jl-1) *  q(jl-1)
                                  * (1.0 + ::mo::constants::c_virtual)
                                  / ((1.0 - q(jl-1) - qcl(jl-1) - qcf(jl-1))
                                      * (1.0 - q(jl-1) - qcl(jl-1)
                                         - qcf(jl-1)))
                                  * vptdrydensHat_jlm1;
          qcfHat(jl-1) += pt(jl-1) *  q(jl-1)
                                  * (1.0 + ::mo::constants::c_virtual)
                                  / ((1.0 - q(jl-1) - qcl(jl-1) - qcf(jl-1))
                                      * (1.0 - q(jl-1) - qcl(jl-1)
                                         - qcf(jl-1)))
                                  * vptdrydensHat_jlm1;
          dryrhoHat(jl) = 0.0;
        }

        exnerHat(0) += dryrho(0) *
          (1.0 - ::mo::constants::rd_over_cp) * dryrhoHat(0) /
          (exner(0) * ::mo::constants::rd_over_cp);

        vptdrydens = pt(0) *
          (1.0 + ::mo::constants::c_virtual * q(0) - qcl(0) - qcf(0)) /
          (1.0 - q(0) - qcl(0) - qcf(0));
        vptdrydensHat = - dryrho(0) * dryrhoHat(0) / vptdrydens;
        ptHat(0) += (1.0 + ::mo::constants::c_virtual * q(0)
                             - qcl(0) - qcf(0))
                            / (1.0 - q(0) - qcl(0) - qcf(0)) * vptdrydensHat;
        qHat(0) += (1.0 + ::mo::constants::c_virtual) * pt(0)
                           * (1.0 - qcl(0) - qcf(0))
                           / ((1.0 - q(0) - qcl(0) - qcf(0))
                           * (1.0 - q(0) - qcl(0) - qcf(0))) * vptdrydensHat;
        qclHat(0) += pt(0) *  q(0) * (1.0 + ::mo::constants::c_virtual)
                             / ((1.0 - q(0) - qcl(0) - qcf(0))
                             * (1.0 - q(0) - qcl(0) - qcf(0))) * vptdrydensHat;
        qcfHat(0) += pt(0) *  q(0) * (1.0 + ::mo::constants::c_virtual)
                             / ((1.0 - q(0) - qcl(0) - qcf(0))
                             * (1.0 - q(0) - qcl(0) - qcf(0))) * vptdrydensHat;
        dryrhoHat(0) = 0.0;
      },
      stateFlds["dry_air_density_levels_minus_one"],
      stateFlds["dimensionless_exner_function_levels_minus_one"],
      stateFlds["height_above_mean_sea_level_levels"],
      stateFlds["height_above_mean_sea_level"],
      stateFlds["air_potential_temperature"],
      stateFlds[specific_humidity_mo],
      stateFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"],
      stateFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"],
      hatFlds["dimensionless_exner_function_levels_minus_one"],
      hatFlds["air_potential_temperature"],
      hatFlds[specific_humidity_mo],
      hatFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"],
      hatFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"],
      hatFlds["dry_air_density_levels_minus_one"]);

  hatFlds["dimensionless_exner_function_levels_minus_one"].set_dirty();
  hatFlds["air_potential_temperature"].set_dirty();
  hatFlds[specific_humidity_mo].set_dirty();
  hatFlds["cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water"].set_dirty();
  hatFlds["cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"].set_dirty();
  hatFlds["dry_air_density_levels_minus_one"].set_dirty();

  oops::Log::trace() << "[eval_dry_air_density_from_exner_levels_minus_one_ad()] ... exit"
                     << std::endl;
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensityFromExnerm1>
  makerDryAirDensityFromExnerm1_("mo_dry_air_density_from_exnerm1");

// -----------------------------------------------------------------------------

DryAirDensityFromExnerm1::DryAirDensityFromExnerm1(const oops::GeometryData & outerGeometryData,
                                                   const oops::Variables & outerVars,
                                                   const eckit::Configuration & covarConf,
                                                   const Parameters_ & params,
                                                   const oops::FieldSet3D & xb,
                                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    xb_(xb.validTime(), xb.commGeom())
{
  xb_.shallowCopy(xb);
  oops::Log::trace() << classname() << "::DryAirDensityFromExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensityFromExnerm1::~DryAirDensityFromExnerm1() {
  oops::Log::trace() << classname() << "::~DryAirDensityFromExnerm1 starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensityFromExnerm1");
  oops::Log::trace() << classname() << "::~DryAirDensityFromExnerm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  eval_dry_air_density_from_exner_levels_minus_one_tl(fset.fieldSet(),
                                                      xb_.fieldSet());
  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_dry_air_density_from_exner_levels_minus_one_ad(fset.fieldSet(),
                                                      xb_.fieldSet());

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  oops::Log::info() << classname() << "::leftInverseMultiply (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  oops::Log::info() << classname() << "::directCalibration (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensityFromExnerm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
