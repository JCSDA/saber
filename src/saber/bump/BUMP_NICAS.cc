/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP_NICAS.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<BUMP_NICAS> makerBUMP_NICAS_("BUMP_NICAS");

// -----------------------------------------------------------------------------

BUMP_NICAS::BUMP_NICAS(const eckit::mpi::Comm & comm,
       const atlas::FunctionSpace & functionSpace,
       const atlas::FieldSet & extraFields,
       const std::vector<size_t> & activeVariableSizes,
       const eckit::Configuration & conf,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg,
       const std::vector<atlas::FieldSet> & fsetVec)
  : SaberCentralBlockBase(conf), bump_()
{
  oops::Log::trace() << classname() << "::BUMP_NICAS starting" << std::endl;

  // Deserialize configuration
  BUMP_NICASParameters params;
  params.validateAndDeserialize(conf);

  // Initialize BUMP
  bump_.reset(new BUMP(comm,
                       functionSpace,
                       extraFields,
                       activeVariableSizes,
                       *params.activeVars.value(),
                       params.bumpParams.value(),
                       fsetVec));

  // Run drivers
  bump_->runDrivers();

  // Partial deallocation
  bump_->partialDealloc();

  oops::Log::trace() << classname() << "::BUMP_NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP_NICAS::~BUMP_NICAS() {
  oops::Log::trace() << classname() << "::~BUMP_NICAS starting" << std::endl;
  util::Timer timer(classname(), "~BUMP_NICAS");
  oops::Log::trace() << classname() << "::~BUMP_NICAS done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_NICAS::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  bump_->randomizeNicas(fset);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_NICAS::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  bump_->multiplyNicas(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP_NICAS::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
