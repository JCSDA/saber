/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP_VerticalBalance.h"

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

static SaberOuterBlockMaker<BUMP_VerticalBalance>
  makerBUMP_VerticalBalance_("BUMP_VerticalBalance");

// -----------------------------------------------------------------------------

BUMP_VerticalBalance::BUMP_VerticalBalance(const eckit::mpi::Comm & comm,
               const oops::GeometryData & outputGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), inputGeometryData_(outputGeometryData), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_VerticalBalance starting"
                     << std::endl;

  // Deserialize configuration
  BUMP_VerticalBalanceParameters params;
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

  oops::Log::trace() << classname() << "::BUMP_VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP_VerticalBalance::~BUMP_VerticalBalance() {
  oops::Log::trace() << classname() << "::~BUMP_VerticalBalance starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_VerticalBalance");
  oops::Log::trace() << classname() << "::~BUMP_VerticalBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_VerticalBalance::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyVbal(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_VerticalBalance::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  bump_->multiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_VerticalBalance::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  bump_->inverseMultiplyVbalAd(fset);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_VerticalBalance::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
