/*
 * (C) Crown Copyright 2022-2024 Met Office
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
#include "mo/eval_total_mixing_ratio.h"
#include "mo/eval_virtual_potential_temperature.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "mo/model2geovals_varchange.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
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
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    xb_(xb.validTime(), xb.commGeom())
{
  xb_.shallowCopy(xb);
  oops::Log::trace() << classname() << "::HydroBal done" << std::endl;
}

// -----------------------------------------------------------------------------

HydroBal::~HydroBal() {
  oops::Log::trace() << classname() << "::~HydroBal starting" << std::endl;
  util::Timer timer(classname(), "~HydroBal");
  oops::Log::trace() << classname() << "::~HydroBal done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_hydrobal_virtual_potential_temperature_tl(fset.fieldSet(), xb_.fieldSet());

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_hydrobal_virtual_potential_temperature_ad(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_hydrobal_hydrostatic_exner_levels_tl(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBal::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
