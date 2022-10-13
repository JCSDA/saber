/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/NICAS.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<NICAS> makerNICAS_("BUMP_NICAS");

// -----------------------------------------------------------------------------

NICAS::NICAS(const oops::GeometryData & geometryData,
             const std::vector<size_t> & activeVariableSizes,
             const oops::Variables & centralVars,
             const Parameters_ & params,
             const atlas::FieldSet & xb,
             const atlas::FieldSet & fg,
             const std::vector<atlas::FieldSet> & fsetVec)
  : bump_()
{
  oops::Log::trace() << classname() << "::NICAS starting" << std::endl;

  // Get active variables
  oops::Variables activeVars = params.activeVars.value().get_value_or(centralVars);

  // Initialize BUMP
  bump_.reset(new BUMP(geometryData.comm(),
                       geometryData.functionSpace(),
                       geometryData.fieldSet(),
                       activeVariableSizes,
                       activeVars,
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

NICAS::~NICAS() {
  oops::Log::trace() << classname() << "::~NICAS starting" << std::endl;
  util::Timer timer(classname(), "~NICAS");
  oops::Log::trace() << classname() << "::~NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  bump_->randomizeNicas(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyNicas(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void NICAS::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
