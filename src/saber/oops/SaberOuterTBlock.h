/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/State.h"
#include "oops/util/FieldSetOperations.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberOuterTBlockParameters : public SaberOuterBlockParametersWrapper {
  OOPS_CONCRETE_PARAMETERS(SaberOuterTBlockParameters, SaberOuterBlockParametersWrapper)
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberOuterTBlock {
  typedef oops::Geometry<MODEL>  Geometry_;
  typedef oops::Increment<MODEL> Increment_;
  typedef oops::State<MODEL>     State_;
  typedef typename std::unordered_map<std::string, const oops::GeometryData*> GeometryDataMap_;

 public:
  typedef SaberOuterTBlockParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::SaberOuterTBlock<MODEL>";}

  SaberOuterTBlock(const Geometry_ &,
                   const GeometryDataMap_ &,
                   const oops::Variables &,
                   const Parameters_ &,
                   const State_ &,
                   const State_ &,
                   const double & adjointTolerance = -1.0);
  ~SaberOuterTBlock() {}

  const GeometryDataMap_ & innerGeometryDataMap() {return innerGeometryDataMap_;}
  const oops::Variables & innerVars() {return saberOuterBlock_->innerVars();}

  void multiply(atlas::FieldSet &) const;
  void multiplyAD(atlas::FieldSet &) const;
  void calibrationInverseMultiply(atlas::FieldSet &) const;

 private:
  void print(std::ostream &);
  GeometryDataMap_ innerGeometryDataMap_;
  std::unique_ptr<SaberOuterBlockBase> saberOuterBlock_;
};

// ----------------------------------------------------------------------------

template <typename MODEL>
SaberOuterTBlock<MODEL>::SaberOuterTBlock(const Geometry_ & geom,
                                          const GeometryDataMap_ & outerGeometryDataMap,
                                          const oops::Variables & outerVars,
                                          const Parameters_ & params,
                                          const State_ & xb,
                                          const State_ & fg,
                                          const double & adjointTolerance)
{
  oops::Log::trace() << classname() << "::SaberOuterTBlock starting" << std::endl;
  // Create outer block
  const SaberBlockParametersBase & saberOuterBlockParams =
    params.saberOuterBlockParameters;
  oops::Log::info() << "Info     : Creating outer block: "
                    << saberOuterBlockParams.saberBlockName.value() << std::endl;

  // Get active variables
  oops::Variables activeVars = saberOuterBlockParams.activeVars.value().get_value_or(outerVars);

  // Read input fields (on model increment geometry)
  std::vector<eckit::LocalConfiguration> inputFields;
  inputFields = saberOuterBlockParams.inputFields.value().get_value_or(inputFields);
  std::vector<atlas::FieldSet> fsetVec = readInputFields(geom,
                                                         activeVars,
                                                         xb.validTime(),
                                                         inputFields);

  // Get active outer variables
  oops::Variables activeOuterVars(outerVars);
  activeOuterVars.intersection(activeVars);

  // Get outer geometryData
  const oops::GeometryData * outerGeometryData = outerGeometryDataMap.at(activeOuterVars[0]);

  // Check that function space type is the same for all active outer variables
  for (const auto var : activeOuterVars.variables()) {
    ASSERT(outerGeometryData->functionSpace().type() ==
      outerGeometryDataMap.at(var)->functionSpace().type());
  }

  // Create outer block
  saberOuterBlock_.reset(SaberOuterBlockFactory::create(
                         *outerGeometryData,
                         geom.variableSizes(activeVars),
                         outerVars,
                         saberOuterBlockParams,
                         xb.fieldSet(),
                         fg.fieldSet(),
                         fsetVec));

  // Inner variables
  const oops::Variables innerVars = saberOuterBlock_->innerVars();

  // Check that active variables are present in either inner or outer variables, or both
  for (const auto & var : activeVars.variables()) {
    ASSERT(innerVars.has(var) || outerVars.has(var));
  }

  // Initialize innerGeometryDataMap
  for (const auto & geometryDataPair : outerGeometryDataMap) {
    innerGeometryDataMap_[geometryDataPair.first] = geometryDataPair.second;
  }

  // Update innerGeometryDataMap with active variables
  for (const auto var : activeVars.variables()) {
    if (innerVars.has(var)) {
      innerGeometryDataMap_[var] = &(saberOuterBlock_->innerGeometryData());
    }
  }

  if (adjointTolerance >= 0.0) {
    // Adjoint test

    // Variables sizes map
    std::unordered_map<std::string, size_t> innerVariableSizeMap;
    for (const auto & var : innerVars.variables()) {
      innerVariableSizeMap[var] = geom.variableSizes(oops::Variables({var}))[0];
    }

    // Create random inner FieldSet
    atlas::FieldSet innerFset = util::createRandomFieldSet(innerGeometryDataMap_,
                                                           innerVariableSizeMap,
                                                           innerVars);

    // Copy inner FieldSet
    atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

    // Variables sizes map
    std::unordered_map<std::string, size_t> outerVariableSizeMap;
    for (const auto & var : outerVars.variables()) {
      outerVariableSizeMap[var] = geom.variableSizes(oops::Variables({var}))[0];
    }

    // Create random outer FieldSet
    atlas::FieldSet outerFset = util::createRandomFieldSet(outerGeometryDataMap,
                                                           outerVariableSizeMap,
                                                           outerVars);

    // Copy outer FieldSet
    atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

    // Apply forward and adjoint multiplication
    saberOuterBlock_->multiply(innerFset);
    saberOuterBlock_->multiplyAD(outerFset);

    // Compute adjoint test
    const double dp1 = util::dotProductFieldSets(innerFset, outerFsetSave, activeVars,
                                                 geom.getComm());
    const double dp2 = util::dotProductFieldSets(outerFset, innerFsetSave, activeVars,
                                                 geom.getComm());
    oops::Log::info() << "Info     : Adjoint test for outer block "
                      << saberOuterBlockParams.saberBlockName.value()
                      << ": y^t (Ax) = " << dp1 << ": x^t (A^t y) = " << dp2 << std::endl;
    ASSERT(abs(dp1) > 0.0);
    ASSERT(abs(dp2) > 0.0);
    oops::Log::test() << "Adjoint test for outer block "
                      << saberOuterBlockParams.saberBlockName.value();
    if (0.5*abs(dp1-dp2)/(dp1+dp2) < adjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      ABORT("Adjoint test failure");
    }
  }

  oops::Log::trace() << classname() << "::SaberOuterTBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterTBlock<MODEL>::multiply(atlas::FieldSet & fset) const
{
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  saberOuterBlock_->multiply(fset);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterTBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const
{
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  saberOuterBlock_->multiplyAD(fset);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterTBlock<MODEL>::calibrationInverseMultiply(atlas::FieldSet & fset) const
{
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;

  saberOuterBlock_->calibrationInverseMultiply(fset);

  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
