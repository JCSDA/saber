/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberOuterBlockBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include <boost/noncopyable.hpp>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberOuterBlockFactory::SaberOuterBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberOuterBlockFactory."
                       << std::endl;
    ABORT("Element already registered in saber::SaberOuterBlockFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

SaberOuterBlockBase * SaberOuterBlockFactory::create(
  const oops::GeometryData & outerGeometryData,
  const std::vector<size_t> & outerVariableSizes,
  const oops::Variables & outerVars,
  const SaberBlockParametersBase & params,
  const atlas::FieldSet & xb,
  const atlas::FieldSet & fg,
  const std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << "SaberOuterBlockBase::create starting" << std::endl;
  const std::string id = params.saberBlockName;
  typename std::map<std::string, SaberOuterBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberOuterBlockFactory." << std::endl;
    ABORT("Element does not exist in saber::SaberOuterBlockFactory.");
  }
  SaberOuterBlockBase * ptr = jsb->second->make(outerGeometryData, outerVariableSizes,
                                                outerVars, params, xb, fg, fsetVec);
  oops::Log::trace() << "SaberOuterBlockBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberBlockParametersBase>
SaberOuterBlockFactory::createParameters(const std::string &name) {
  typename std::map<std::string, SaberOuterBlockFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in saber::SaberOuterBlockFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

void SaberOuterBlockBase::adjointTest(const eckit::mpi::Comm & comm,
                                      const oops::GeometryData & outerGeometryData,
                                      const std::vector<size_t> & outerVariableSizes,
                                      const oops::Variables & outerVars,
                                      const oops::GeometryData & innerGeometryData,
                                      const std::vector<size_t> & innerVariableSizes,
                                      const oops::Variables & innerVars,
                                      const double & adjointTolerance) const {
  oops::Log::trace() << "SaberOuterBlockBase::adjointTest starting" << std::endl;

  // Create random inner FieldSet
  atlas::FieldSet innerFset = util::createRandomFieldSet(innerGeometryData,
                                                         innerVariableSizes,
                                                         innerVars);

  // Copy inner FieldSet
  atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

  // Create random outer FieldSet
  atlas::FieldSet outerFset = util::createRandomFieldSet(outerGeometryData,
                                                         outerVariableSizes,
                                                         outerVars);

  // Copy outer FieldSet
  atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

  // Apply forward and adjoint multiplication
  this->multiply(innerFset);
  this->multiplyAD(outerFset);

  // Compute adjoint test
  const double dp1 = util::dotProductFieldSets(innerFset, outerFsetSave, outerVars, comm);
  const double dp2 = util::dotProductFieldSets(outerFset, innerFsetSave, innerVars, comm);
  oops::Log::info() << "Info     : Adjoint test: y^t (Ax) = " << dp1
                    << ": x^t (A^t y) = " << dp2 << std::endl;
  oops::Log::test() << "Adjoint test";
  if (0.5*abs(dp1-dp2)/(dp1+dp2) < adjointTolerance) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    ABORT("Adjoint test failure");
  }

  oops::Log::trace() << "SaberOuterBlockBase::adjointTest done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
