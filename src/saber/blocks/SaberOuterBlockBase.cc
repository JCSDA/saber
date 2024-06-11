/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberOuterBlockBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
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

#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberOuterBlockFactory::SaberOuterBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberOuterBlockFactory."
                       << std::endl;
    throw eckit::Exception("Element already registered in saber::SaberOuterBlockFactory.", Here());
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberOuterBlockBase> SaberOuterBlockFactory::create(
  const oops::GeometryData & outerGeometryData,
  const oops::Variables & outerVars,
  const eckit::Configuration & covarConfig,
  const SaberBlockParametersBase & params,
  const oops::FieldSet3D & xb,
  const oops::FieldSet3D & fg) {
  oops::Log::trace() << "SaberOuterBlockBase::create starting" << std::endl;
  const std::string id = params.saberBlockName;
  typename std::map<std::string, SaberOuterBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberOuterBlockFactory." << std::endl;
    throw eckit::UserError("Element does not exist in saber::SaberOuterBlockFactory.", Here());
  }
  std::unique_ptr<SaberOuterBlockBase> ptr =
    jsb->second->make(outerGeometryData, outerVars, covarConfig, params, xb, fg);
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

void SaberOuterBlockBase::adjointTest(const oops::GeometryData & outerGeometryData,
                                      const oops::Variables & outerVars,
                                      const oops::GeometryData & innerGeometryData,
                                      const oops::Variables & innerVars,
                                      const double & adjointTolerance) const {
  oops::Log::trace() << "SaberOuterBlockBase::adjointTest starting" << std::endl;

  // Create random inner FieldSet
  oops::FieldSet3D innerFset = oops::randomFieldSet3D(validTime_,
                                                      innerGeometryData.comm(),
                                                      innerGeometryData.functionSpace(),
                                                      innerVars);

  // Copy inner FieldSet
  oops::FieldSet3D innerFsetSave(innerFset);

  // Create random outer FieldSet
  oops::FieldSet3D outerFset = oops::randomFieldSet3D(validTime_,
                                                      outerGeometryData.comm(),
                                                      outerGeometryData.functionSpace(),
                                                      outerVars);

  // Set its halos to zero, as they should be coming in to the block
  if (outerGeometryData.functionSpace().type() != "Spectral") {
    auto ghostView = atlas::array::make_view<int, 1>(outerGeometryData.functionSpace().ghost());
    for (const auto & var : outerVars) {
      auto view = atlas::array::make_view<double, 2>(outerFset[var]);
      for (int jn = 0; jn < view.shape(0); ++jn) {
        if (ghostView(jn) == 1) {
          for (int jl = 0; jl < view.shape(1); ++jl) {
            view(jn, jl) = 0.0;
          }
        }
      }
    }
  }

  // Copy outer FieldSet
  oops::FieldSet3D outerFsetSave(outerFset);

  // Apply forward and adjoint multiplication
  this->multiply(innerFset);
  this->multiplyAD(outerFset);

  // Compute adjoint test
  const double dp1 = innerFset.dot_product_with(outerFsetSave, outerVars);
  const double dp2 = outerFset.dot_product_with(innerFsetSave, innerVars);
  oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                    << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                    << adjointTolerance << std::endl;
  oops::Log::test() << "Adjoint test for block " << this->blockName();
  if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < adjointTolerance) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    throw eckit::Exception("Adjoint test failure for block " + this->blockName(), Here());
  }

  oops::Log::trace() << "SaberOuterBlockBase::adjointTest done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberOuterBlockBase::inverseTest(const oops::GeometryData & innerGeometryData,
                                      const oops::Variables & innerVars,
                                      const oops::GeometryData & outerGeometryData,
                                      const oops::Variables & outerVars,
                                      const oops::Variables & innerVarsToCompare,
                                      const oops::Variables & outerVarsToCompare,
                                      const double & innerInverseTolerance,
                                      const double & outerInverseTolerance) const {
  oops::Log::trace() << "SaberOuterBlockBase::inverseTest starting" << std::endl;

  // Inner inverse test

  // Create inner FieldSet
  oops::FieldSet3D innerFset = this->generateInnerFieldSet(innerGeometryData,
                                                           innerVars);

  // Apply forward multiplication
  this->multiply(innerFset);

  // Save inner FieldSet
  oops::FieldSet3D innerFsetSave(innerFset);

  // Apply inverse multiplication
  this->leftInverseMultiply(innerFset);

  // Apply forward multiplication
  this->multiply(innerFset);

  // Check that the fieldsets contain the same fields
  auto innerFieldNames = oops::Variables(innerFset.field_names());
  auto innerFieldNamesSave = oops::Variables(innerFsetSave.field_names());

  if (innerFieldNames != innerFieldNamesSave) {
    throw eckit::Exception("Inner inverse test for block " + this->blockName()
      + ": fieldsets content does not match", Here());
  }

  // Check that the fieldsets are similar within tolerance
  oops::Variables outerVariablesToRemove(innerFieldNames);
  outerVariablesToRemove -= outerVarsToCompare;

  innerFset.removeFields(outerVariablesToRemove);
  innerFsetSave.removeFields(outerVariablesToRemove);
  const bool outerComparison = this->compareFieldSets(innerFset,
                                                      innerFsetSave,
                                                      innerInverseTolerance);
  oops::Log::test() << "Inner inverse test for block " << this->blockName();
  if (outerComparison) {
    oops::Log::test() << " passed: U Uinv (U x) == (U x)" << std::endl;
  } else {
    oops::Log::test() << " failed: U Uinv (U x) != (U x)" << std::endl;
    throw eckit::Exception("Inner inverse test failure for block " + this->blockName(), Here());
  }

  // Outer inverse test

  // Create outer FieldSet
  oops::FieldSet3D outerFset = this->generateOuterFieldSet(outerGeometryData,
                                                           outerVars);

  // Apply inverse multiplication
  this->leftInverseMultiply(outerFset);

  // Save outer FieldSet
  oops::FieldSet3D outerFsetSave(outerFset);

  // Apply forward multiplication
  this->multiply(outerFset);

  // Apply inverse multiplication
  this->leftInverseMultiply(outerFset);

  // Check that the fieldsets contain the same fields
  auto outerFieldNames = oops::Variables(outerFset.field_names());
  auto outerFieldNamesSave = oops::Variables(outerFsetSave.field_names());
  if (outerFieldNames != outerFieldNamesSave) {
    throw eckit::Exception("Outer inverse test for block " + this->blockName()
      + ": fieldsets content does not match", Here());
  }

  // Check that the fieldsets are similar within tolerance
  oops::Variables innerVariablesToRemove(outerFieldNames);
  innerVariablesToRemove -= innerVarsToCompare;
  outerFset.removeFields(innerVariablesToRemove);
  outerFsetSave.removeFields(innerVariablesToRemove);
  const bool innerComparison = this->compareFieldSets(outerFset,
                                                      outerFsetSave,
                                                      outerInverseTolerance);
  oops::Log::test() << "Outer inverse test for block " << this->blockName();
  if (innerComparison) {
    oops::Log::test() << " passed: Uinv U (Uinv x) == (Uinv x)" << std::endl;
  } else {
    oops::Log::test() << " failed: Uinv U (Uinv x) != (Uinv x)" << std::endl;
    throw eckit::Exception("Outer inverse test failure for block " + this->blockName(), Here());
  }

  oops::Log::trace() << "SaberOuterBlockBase::inverseTest done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
