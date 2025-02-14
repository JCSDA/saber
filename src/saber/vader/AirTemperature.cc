/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/AirTemperature.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

// note below routine this is also used in H(x), so we keep it in vader
#include "mo/eval_air_temperature.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<AirTemperature> makerAirTemperature_("mo_air_temperature");

// -----------------------------------------------------------------------------

AirTemperature::AirTemperature(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::AirTemperature done" << std::endl;
}

// -----------------------------------------------------------------------------

AirTemperature::~AirTemperature() {
  oops::Log::trace() << classname() << "::~AirTemperature starting" << std::endl;
  util::Timer timer(classname(), "~AirTemperature");
  oops::Log::trace() << classname() << "::~AirTemperature done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperature::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  mo::eval_air_temperature_tl(fset.fieldSet(), xb_.fieldSet());

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperature::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  mo::eval_air_temperature_ad(fset.fieldSet(), xb_.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void AirTemperature::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
