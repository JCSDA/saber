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
#include "oops/util/FieldSetHelpers.h"
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
  const eckit::Configuration & covarConfig,
  const SaberBlockParametersBase & params,
  const atlas::FieldSet & xb,
  const atlas::FieldSet & fg,
  const util::DateTime & validTime) {
  oops::Log::trace() << "SaberOuterBlockBase::create starting" << std::endl;
  const std::string id = params.saberBlockName;
  typename std::map<std::string, SaberOuterBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberOuterBlockFactory." << std::endl;
    ABORT("Element does not exist in saber::SaberOuterBlockFactory.");
  }
  SaberOuterBlockBase * ptr = jsb->second->make(outerGeometryData, outerVariableSizes,
                                                outerVars, covarConfig, params, xb, fg, validTime);
  ptr->setBlockName(params.saberBlockName.value());
  ptr->setSkipInverse(params.skipInverse.value());
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
  atlas::FieldSet innerFset = util::createRandomFieldSet(innerGeometryData.comm(),
                                                         innerGeometryData.functionSpace(),
                                                         innerVariableSizes,
                                                         innerVars.variables());

  // Copy inner FieldSet
  atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

  // Create random outer FieldSet
  atlas::FieldSet outerFset = util::createRandomFieldSet(outerGeometryData.comm(),
                                                         outerGeometryData.functionSpace(),
                                                         outerVariableSizes,
                                                         outerVars.variables());

  // Copy outer FieldSet
  atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

  // Apply forward and adjoint multiplication
  this->multiply(innerFset);
  this->multiplyAD(outerFset);

  // Compute adjoint test
  const double dp1 = util::dotProductFieldSets(innerFset, outerFsetSave,
                                               outerVars.variables(), comm);
  const double dp2 = util::dotProductFieldSets(outerFset, innerFsetSave,
                                               innerVars.variables(), comm);
  oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                    << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                    << adjointTolerance << std::endl;
  oops::Log::test() << "Adjoint test for block " << this->blockName();
  if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < adjointTolerance) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    ABORT("Adjoint test failure for block " + this->blockName());
  }

  oops::Log::trace() << "SaberOuterBlockBase::adjointTest done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberOuterBlockBase::inverseTest(const oops::GeometryData & innerGeometryData,
                                      const std::vector<size_t> & innerVariableSizes,
                                      const oops::Variables & innerVars,
                                      const oops::GeometryData & outerGeometryData,
                                      const std::vector<size_t> & outerVariableSizes,
                                      const oops::Variables & outerVars,
                                      const oops::Variables & innerVarsToCompare,
                                      const oops::Variables & outerVarsToCompare,
                                      const double & innerInverseTolerance,
                                      const double & outerInverseTolerance,
                                      const size_t & timeRank) const {
  oops::Log::trace() << "SaberOuterBlockBase::inverseTest starting" << std::endl;

  // Inner inverse test

  // Create inner FieldSet
  atlas::FieldSet innerFset = this->generateInnerFieldSet(innerGeometryData,
                                                          innerVariableSizes,
                                                          innerVars,
                                                          timeRank);

  // Apply forward multiplication
  this->multiply(innerFset);

  // Save inner FieldSet
  atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

  // Apply inverse multiplication
  this->leftInverseMultiply(innerFset);

  // Apply forward multiplication
  this->multiply(innerFset);

  // Check that the fieldsets contain the same fields
  auto innerFieldNames = oops::Variables(innerFset.field_names());
  auto innerFieldNamesSave = oops::Variables(innerFsetSave.field_names());

  if (innerFieldNames != innerFieldNamesSave) {
    ABORT("Inner inverse test for block " + this->blockName()
      + ": fieldsets content does not match");
  }

  // Check that the fieldsets are similar within tolerance
  oops::Variables outerVariablesToRemove(innerFieldNames);
  outerVariablesToRemove -= outerVarsToCompare;

  util::removeFieldsFromFieldSet(innerFset, outerVariablesToRemove.variables());
  util::removeFieldsFromFieldSet(innerFsetSave, outerVariablesToRemove.variables());
  const bool outerComparison = this->compareFieldSets(innerFset,
                                                      innerFsetSave,
                                                      innerInverseTolerance);
  oops::Log::test() << "Inner inverse test passed for block " << this->blockName();
  if (outerComparison) {
    oops::Log::test() << " passed: U Uinv (U x) == (U x)" << std::endl;
  } else {
    oops::Log::test() << " failed: U Uinv (U x) != (U x)" << std::endl;
    ABORT("Inner inverse test failure for block " + this->blockName());
  }

  // Outer inverse test

  // Create outer FieldSet
  atlas::FieldSet outerFset = this->generateOuterFieldSet(outerGeometryData,
                                                          outerVariableSizes,
                                                          outerVars,
                                                          timeRank);

  // Apply inverse multiplication
  this->leftInverseMultiply(outerFset);

  // Save outer FieldSet
  atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

  // Apply forward multiplication
  this->multiply(outerFset);

  // Apply inverse multiplication
  this->leftInverseMultiply(outerFset);

  // Check that the fieldsets contain the same fields
  auto outerFieldNames = oops::Variables(outerFset.field_names());
  auto outerFieldNamesSave = oops::Variables(outerFsetSave.field_names());
  if (outerFieldNames != outerFieldNamesSave) {
    ABORT("Outer inverse test for block " + this->blockName()
      + ": fieldsets content does not match");
  }

  // Check that the fieldsets are similar within tolerance
  oops::Variables innerVariablesToRemove(outerFieldNames);
  innerVariablesToRemove -= innerVarsToCompare;
  util::removeFieldsFromFieldSet(outerFset, innerVariablesToRemove.variables());
  util::removeFieldsFromFieldSet(outerFsetSave, innerVariablesToRemove.variables());
  const bool innerComparison = this->compareFieldSets(outerFset,
                                                      outerFsetSave,
                                                      outerInverseTolerance);
  oops::Log::test() << "Outer inverse test for block " << this->blockName();
  if (innerComparison) {
    oops::Log::test() << " passed: Uinv U (Uinv x) == (Uinv x)" << std::endl;
  } else {
    oops::Log::test() << " failed: Uinv U (Uinv x) != (Uinv x)" << std::endl;
    ABORT("Outer inverse test failure for block " + this->blockName());
  }

  oops::Log::trace() << "SaberOuterBlockBase::inverseTest done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
