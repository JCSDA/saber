/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/StdDev.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<StdDev> makerStdDev_("BUMP_StdDev");

// -----------------------------------------------------------------------------

StdDev::StdDev(const oops::GeometryData & outputGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const oops::Variables & outputVars,
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : inputGeometryData_(outputGeometryData), inputVars_(outputVars), bump_()
{
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

  // Get active variables
  oops::Variables activeVars = params.activeVars.value().get_value_or(outputVars);

  // Initialize BUMP
  bump_.reset(new BUMP(outputGeometryData.comm(),
                       outputGeometryData.functionSpace(),
                       outputGeometryData.fieldSet(),
                       activeVariableSizes,
                       activeVars,
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

StdDev::~StdDev() {
  oops::Log::trace() << classname() << "::~StdDev starting" << std::endl;
  util::Timer timer(classname(), "~StdDev");
  oops::Log::trace() << classname() << "::~StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyStdDev(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyStdDev(fset);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
