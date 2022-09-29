/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP_StdDev.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BUMP_StdDev> makerBUMP_StdDev_("BUMP_StdDev");

// -----------------------------------------------------------------------------

BUMP_StdDev::BUMP_StdDev(const eckit::mpi::Comm & comm,
               const oops::GeometryData & outputGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), inputGeometryData_(outputGeometryData), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_StdDev starting" << std::endl;

  // Deserialize configuration
  BUMP_StdDevParameters params;
  params.deserialize(conf);

  // Input variables
  inputVars_ = params.outputVars.value();

  // Initialize BUMP
  bump_.reset(new BUMP(comm,
                       outputGeometryData.functionSpace(),
                       outputGeometryData.fieldSet(),
                       activeVariableSizes,
                       *params.activeVars.value(),
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::BUMP_StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP_StdDev::~BUMP_StdDev() {
  oops::Log::trace() << classname() << "::~BUMP_StdDev starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_StdDev");
  oops::Log::trace() << classname() << "::~BUMP_StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyStdDev(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyStdDev(fset);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
