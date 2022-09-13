/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/interpolation/GSI_Interpolation.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/interpolation/GSI_InterpolationImpl.h"
#include "saber/oops/SaberBlockBase.h"


namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberBlockMaker<gsi::Interpolation>
  makerGSI_Interpolation_("gsi interpolation to model grid");

// -------------------------------------------------------------------------------------------------

Interpolation::Interpolation(const eckit::mpi::Comm & comm,
                             const atlas::FunctionSpace & functionSpace,
                             const atlas::FieldSet & extraFields,
                             const std::vector<size_t> & variableSizes,
                             const Parameters_ & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params),
interpolationImpl_(comm, functionSpace, params, params.inputVars.value().variables())
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;
  util::Timer timer(classname(), "Interpolation");
  // Assert that there is no variable change in this block
  ASSERT(params.inputVars.value() == params.outputVars.value());
  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Interpolation::~Interpolation() {
  oops::Log::trace() << classname() << "::~Interpolation starting" << std::endl;
  util::Timer timer(classname(), "~Interpolation");
  oops::Log::trace() << classname() << "::~Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::randomize(atlas::FieldSet & fset) const {
  ABORT(classname() + "randomize: not implemented");
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  interpolationImpl_.multiply(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::inverseMultiply(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiply: not implemented");
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  interpolationImpl_.multiplyAD(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::inverseMultiplyAD(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiplyAD: not implemented");
}

// -------------------------------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
