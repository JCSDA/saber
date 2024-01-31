/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HpHexnerToPExner.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HpHexnerToPExner>
  makerHpHexnerToPExner_(
  "mo_hydrostatic_pressure_hydrostatic_exner_to_pressure_exner");

// -----------------------------------------------------------------------------

HpHexnerToPExner::HpHexnerToPExner(const oops::GeometryData & outerGeometryData,
                     const oops::Variables & outerVars,
                     const eckit::Configuration & covarConf,
                     const Parameters_ & params,
                     const oops::FieldSet3D & xb,
                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars))
{
  oops::Log::trace() << classname() << "::HpHexnerToPExner starting" << std::endl;
  oops::Log::trace() << classname() << "::HpHexnerToPExner done" << std::endl;
}

// -----------------------------------------------------------------------------

HpHexnerToPExner::~HpHexnerToPExner() {
  oops::Log::trace() << classname() << "::~HpHexnerToPExner starting" << std::endl;
  util::Timer timer(classname(), "~HpHexnerToPExner");
  oops::Log::trace() << classname() << "::~HpHexnerToPExner done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpHexnerToPExner::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  auto hydrostaticPressureView =
    atlas::array::make_view<const double, 2>(fset["hydrostatic_pressure_levels"]);
  auto hydrostaticExnerView =
    atlas::array::make_view<const double, 2>(fset["hydrostatic_exner_levels"]);
  auto airPressureView =
    atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  auto exnerLevelsMinusOneView =
     atlas::array::make_view<double, 2>(fset["exner_levels_minus_one"]);

  airPressureView.assign(hydrostaticPressureView);
  // Note that the number of levels in hydrostaticExnerView is one more than
  // in exnerLevelsMinusOneView. This is o.k. however as
  // the assign method copies the region that is common to both Views.
  exnerLevelsMinusOneView.assign(hydrostaticExnerView);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpHexnerToPExner::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  auto airPressureView =
      atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  auto hydrostaticPressureView =
      atlas::array::make_view<double, 2>(fset["hydrostatic_pressure_levels"]);
  auto hydrostaticExnerView =
      atlas::array::make_view<double, 2>(fset["hydrostatic_exner_levels"]);
  auto exnerLevelsMinusOneView =
      atlas::array::make_view<double, 2>(fset["exner_levels_minus_one"]);

  for (atlas::idx_t jn = 0; jn < fset["exner_levels_minus_one"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < fset["exner_levels_minus_one"].shape(1); ++jl) {
      hydrostaticExnerView(jn, jl) += exnerLevelsMinusOneView(jn, jl);
      exnerLevelsMinusOneView(jn, jl) = 0.0;
    }
  }
  for (atlas::idx_t jn = 0; jn < fset["hydrostatic_pressure_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < fset["hydrostatic_pressure_levels"].shape(1); ++jl) {
      hydrostaticPressureView(jn, jl) += airPressureView(jn, jl);
      airPressureView(jn, jl) = 0.0;
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpHexnerToPExner::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  // Retrieve hydrostatic Exner from Exner. Need to extrapolate top level
  auto exner_view = atlas::array::make_view<const double, 2>(fset["exner_levels_minus_one"]);
  auto hexner_view = atlas::array::make_view<double, 2>(fset["hydrostatic_exner_levels"]);
  auto airPressureView =
    atlas::array::make_view<const double, 2>(fset["air_pressure_levels"]);
  auto hydrostaticPressureView =
    atlas::array::make_view<double, 2>(fset["hydrostatic_pressure_levels"]);

  hydrostaticPressureView.assign(airPressureView);

  const auto levels = fset["hydrostatic_exner_levels"].shape(1);
  for (atlas::idx_t jnode = 0; jnode < hexner_view.shape(0); jnode++) {
    for (atlas::idx_t jlev = 0; jlev < levels - 1; jlev++) {
      hexner_view(jnode, jlev) = exner_view(jnode, jlev);
    }
    // Extrapolate to top level assuming the same hydrostatic Exner increment as in level below.
    // This is consistent with the extrapolation of hydrostatic pressure increment done in
    // mo::eval_hydrostatic_exner_tl.
    hexner_view(jnode, levels - 1) = hexner_view(jnode, levels - 2);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HpHexnerToPExner::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
