/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/VerticalBalance.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VerticalBalance> makerVerticalBalance_("BUMP_VerticalBalance");

// -----------------------------------------------------------------------------

VerticalBalance::VerticalBalance(const oops::GeometryData & outerGeometryData,
                                 const std::vector<size_t> & activeVariableSizes,
                                 const oops::Variables & outerVars,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const std::vector<atlas::FieldSet> & fsetVec)
  : innerGeometryData_(outerGeometryData), innerVars_(outerVars), bump_()
{
  oops::Log::trace() << classname() << "::VerticalBalance starting"
                     << std::endl;

  // Get active variables
  oops::Variables activeVars = params.activeVars.value().get_value_or(outerVars);

  // Check active variables sizes vector consistency
  if (std::adjacent_find(activeVariableSizes.begin(), activeVariableSizes.end(),
    std::not_equal_to<>()) != activeVariableSizes.end()) {
    ABORT("inconsistent elements in activeVariableSizes");
  }

  // Initialize BUMP
  bump_.reset(new BUMP(outerGeometryData.comm(),
                       outerGeometryData.functionSpace(),
                       outerGeometryData.fieldSet(),
                       activeVariableSizes,
                       activeVars,
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

VerticalBalance::~VerticalBalance() {
  oops::Log::trace() << classname() << "::~VerticalBalance starting" << std::endl;
  util::Timer timer(classname(), "~VerticalBalance");
  oops::Log::trace() << classname() << "::~VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyVbal(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyVbal(fset);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VerticalBalance::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
