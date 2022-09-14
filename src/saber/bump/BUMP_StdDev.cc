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
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<BUMP_StdDev> makerBUMP_StdDev_("BUMP_StdDev");

// -----------------------------------------------------------------------------

BUMP_StdDev::BUMP_StdDev(const eckit::mpi::Comm & comm,
                         const atlas::FunctionSpace & functionSpace,
                         const atlas::FieldSet & extraFields,
                         const std::vector<size_t> & variableSizes,
                         const Parameters_ & params,
                         const atlas::FieldSet & xb,
                         const atlas::FieldSet & fg,
                         const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_StdDev starting" << std::endl;

  // Setup and check input/ouput variables
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  ASSERT(inputVars == outputVars);

  // Active variables
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= inputVars);
  } else {
    activeVars += inputVars;
  }

  // Initialize BUMP
  bump_.reset(new BUMP(comm,
                       functionSpace,
                       extraFields,
                       variableSizes,
                       activeVars,
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

void BUMP_StdDev::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyStdDev(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  bump_->inverseMultiplyStdDev(fset);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  this->inverseMultiply(fset);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
