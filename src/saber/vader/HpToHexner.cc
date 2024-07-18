/*
 * (C) Crown Copyright 2022-2024 Met Office
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

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HpToHexner>
  makerHpToHexner_("mo_hydrostatic_pressure_to_hydrostatic_exner");

// -----------------------------------------------------------------------------

HpToHexner::HpToHexner(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(getActiveVars(params, outerVars)),
    xb_(xb)
{
  oops::Log::trace() << classname() << "::HpToHexner done" << std::endl;
}

// -----------------------------------------------------------------------------

HpToHexner::~HpToHexner() {
  oops::Log::trace() << classname() << "::~HpToHexner starting" << std::endl;
  util::Timer timer(classname(), "~HpToHexner");
  oops::Log::trace() << classname() << "::~HpToHexner done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  const oops::Variables outputVars({"hydrostatic_exner_levels"});
  allocateMissingFields(fset,
                        outputVars,
                        activeVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_hydrostatic_exner_levels_tl(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::eval_hydrostatic_exner_levels_ad(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Retrieve hydrostatic pressure from hydrostatic Exner.
  mo::eval_hydrostatic_exner_levels_tl_inv(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpToHexner::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
