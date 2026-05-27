/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryMoistIncrOp.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_cloud_liquid_water_mixing_ratio.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "oops/base/Variables.h"
#include "oops/util/for_each.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryMoistIncrOp>
  makerDryMoistIncrOp_("mo_dry_mio");

// -----------------------------------------------------------------------------

DryMoistIncrOp::DryMoistIncrOp(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::DryMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

DryMoistIncrOp::~DryMoistIncrOp() {
  oops::Log::trace() << classname() << "::~DryMoistIncrOp starting" << std::endl;
  util::Timer timer(classname(), "~DryMoistIncrOp");
  oops::Log::trace() << classname() << "::~DryMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryMoistIncrOp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());
  // m_v
  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_inv_tl(
    fset.fieldSet(), xb_.fieldSet());

  // m_cl
  mo::eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_inv_tl(
    fset.fieldSet(), xb_.fieldSet());

  // m_ci
  mo::eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_inv_tl(
    fset.fieldSet(), xb_.fieldSet());

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryMoistIncrOp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_inv_ad(
    fset.fieldSet(), xb_.fieldSet());

  mo::eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_inv_ad(
    fset.fieldSet(), xb_.fieldSet());

  mo::eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_inv_ad(
    fset.fieldSet(), xb_.fieldSet());

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryMoistIncrOp::inverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
    fset.fieldSet(), xb_.fieldSet());

  mo::eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
    fset.fieldSet(), xb_.fieldSet());

  mo::eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_tl(
    fset.fieldSet(), xb_.fieldSet());

  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryMoistIncrOp::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  oops::Log::info() << classname() << "::directCalibration (empty step)" << std::endl;
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------
void DryMoistIncrOp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
