/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/DryAirDensity.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_cloud_liquid_mixing_ratio.h"
#include "mo/eval_dry_air_density.h"
#include "mo/eval_total_mixing_ratio.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

using atlas::array::make_view;
using atlas::idx_t;

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DryAirDensity> makerDryAirDensity_("mo_dry_air_density");

// -----------------------------------------------------------------------------

DryAirDensity::DryAirDensity(const oops::GeometryData & outerGeometryData,
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
    intermediateTempVars_(params.intermediateTempVars(outerVars)),
    xb_(xb.validTime(), xb.commGeom())
{
  xb_.shallowCopy(xb);
  oops::Log::trace() << classname() << "::DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

DryAirDensity::~DryAirDensity() {
  oops::Log::trace() << classname() << "::~DryAirDensity starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensity");
  oops::Log::trace() << classname() << "::~DryAirDensity done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());
  // Allocate temporary pressure_levels_minus_one field
  allocateMissingFields(fset, intermediateTempVars_, intermediateTempVars_,
                        innerGeometryData_.functionSpace());

  // Populate pressure_levels_minus_one from pressure_levels, assumed already populated
  const auto pIncView = make_view<const double, 2>(fset["air_pressure_levels"]);
  auto pmoIncView = make_view<double, 2>(fset["air_pressure_levels_minus_one"]);

  pmoIncView.assign(pIncView);

  // Populate output fields.
  mo::eval_dry_air_density_from_pressure_levels_minus_one_tl(fset.fieldSet(),
                                                             xb_.fieldSet());

  // Remove temporary pressure_levels_minus_one field
  fset.removeFields(intermediateTempVars_);
  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------

void DryAirDensity::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());
  // Allocate temporary pressure_levels_minus_one field
  allocateMissingFields(fset, intermediateTempVars_, intermediateTempVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_dry_air_density_from_pressure_levels_minus_one_ad(fset.fieldSet(),
                                                             xb_.fieldSet());

  // Populate pressure_levels from pressure_levels_minus_one.
  // Note there is no update to pHatView(jn,top_jl),
  // where top_jl = fset["air_pressure_levels_minus_one"].shape(1)
  // is the uppermost air_pressure_levels value,
  // as eval_dry_air_density_from_pressure_levels_minus_one_ad updates the adjoint
  // of dry_air_density_from_pressure_levels_minus_one, which does not have this
  // uppermost level
  auto pHatView = make_view<double, 2>(fset["air_pressure_levels"]);
  auto pmoHatView = make_view<double, 2>(fset["air_pressure_levels_minus_one"]);

  for (idx_t jn = 0; jn < fset["air_pressure_levels_minus_one"].shape(0); ++jn) {
    for (idx_t jl = 0; jl < fset["air_pressure_levels_minus_one"].shape(1); ++jl) {
      pHatView(jn, jl) += pmoHatView(jn, jl);
      pmoHatView(jn, jl) = 0.0;
    }
  }

  // Remove temporary pressure_levels_minus_one field
  fset.removeFields(intermediateTempVars_);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void DryAirDensity::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
