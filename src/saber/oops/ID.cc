/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/ID.h"

#include <memory>
#include <string>
#include <vector>

#include "oops/util/Timer.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<ID> makerID_("ID");

// -----------------------------------------------------------------------------

ID::ID(const atlas::FunctionSpace & functionSpace,
       const atlas::FieldSet & extraFields,
       const std::vector<size_t> & variableSizes,
       const Parameters_ & params,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg,
       const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params)
{
  oops::Log::trace() << classname() << "::ID starting" << std::endl;
  oops::Log::trace() << classname() << "::ID done" << std::endl;
}

// -----------------------------------------------------------------------------

ID::~ID() {
  oops::Log::trace() << classname() << "::~ID starting" << std::endl;
  util::Timer timer(classname(), "~ID");
  oops::Log::trace() << classname() << "::~ID done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
