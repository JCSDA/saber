/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/StdDev.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<StdDev> makerStdDev_("StdDev");

// -----------------------------------------------------------------------------

StdDev::StdDev(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & inputFunctionSpace,
               const atlas::FieldSet & inputExtraFields,
               const std::vector<size_t> & inputVariableSizes,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & outputVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), stdDevFset_()
{
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

  // Deserialize configuration
  StdDevParameters params;
  params.deserialize(conf);

  // Copy stddev field
  stdDevFset_.clear();
  for (const auto & fset : fsetVec) {
    if (fset.name() == "StdDev") {
      for (const auto & field : fset) { // TODO: should it be only active variables?
        stdDevFset_.add(field);
      }
    }
  }

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
  util::FieldSetMultiply(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::FieldSetMultiply(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  util::FieldSetDivide(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
