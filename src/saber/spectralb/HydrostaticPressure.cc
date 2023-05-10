/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/HydrostaticPressure.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/spectralb/GaussUVToGP.h"
#include "saber/vader/GpToHp.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydrostaticPressure>
  makerHydrostaticPressure_("mo_hydrostatic_pressure");

// -----------------------------------------------------------------------------

HydrostaticPressure::HydrostaticPressure(const oops::GeometryData & outerGeometryData,
                                   const std::vector<size_t> & activeVariableSizes,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const atlas::FieldSet & xb,
                                   const atlas::FieldSet & fg,
                                   const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    gptohp_(std::make_unique<saber::vader::GpToHp>(outerGeometryData,
                                                   activeVariableSizes,
                                                   outerVars,
                                                   covarConf,
                                                   params.gpToHp,
                                                   xb, fg, validTimeOfXbFg)),
    gaussuvtogp_(std::make_unique<GaussUVToGP>(outerGeometryData,
                                               activeVariableSizes,
                                               gptohp_->innerVars(),
                                               covarConf,
                                               params.gaussUVToGp, xb, fg, validTimeOfXbFg))
{
  oops::Log::trace() << classname() << "::HydrostaticPressure starting" << std::endl;
  oops::Log::trace() << classname() << "::HydrostaticPressure done" << std::endl;
}

// -----------------------------------------------------------------------------

HydrostaticPressure::~HydrostaticPressure() {
  oops::Log::trace() << classname() << "::~HydrostaticPressure starting" << std::endl;
  util::Timer timer(classname(), "~HydrostaticPressure");
  oops::Log::trace() << classname() << "::~HydrostaticPressure done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticPressure::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // note gaussuvtogp_->multiply creates "geostrophic_pressure_levels_minus_one"
  // if not there.
  gaussuvtogp_->multiply(fset);
  gptohp_->multiply(fset);
  // remove "geostrophic_pressure_levels_minus_one" since it is not an
  // active variable (but a temporary one).
  util::removeFieldsFromFieldSet(fset, {"geostrophic_pressure_levels_minus_one"});
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticPressure::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // add "geostrophic_pressure_levels_minus_one" to fset and assign to zero

  atlas::Field gp = gaussFunctionSpace_.createField<double>(
    atlas::option::name("geostrophic_pressure_levels_minus_one") |
    atlas::option::levels(fset["eastward_wind"].levels()));
  gp.haloExchange();
  atlas::array::make_view<double, 2>(gp).assign(0.0);
  fset.add(gp);

  gptohp_->multiplyAD(fset);
  gaussuvtogp_->multiplyAD(fset);
  // remove "geostrophic_pressure_levels_minus_one" since it is not an
  // active variable (but a temporary one)."
  util::removeFieldsFromFieldSet(fset, {"geostrophic_pressure_levels_minus_one"});
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticPressure::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  gaussuvtogp_->multiply(fset);
  gptohp_->leftInverseMultiply(fset);
  util::removeFieldsFromFieldSet(fset, {"geostrophic_pressure_levels_minus_one"});
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticPressure::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
