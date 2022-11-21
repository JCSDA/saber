/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
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
  typedef oops::Geometry<MODEL>                              Geometry_;
  typedef oops::Increment<MODEL>                             Increment_;
  typedef oops::State<MODEL>                                 State_;

 public:
  typedef SaberOuterTBlockParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::SaberOuterTBlock<MODEL>";}

  SaberOuterTBlock(const Geometry_ &,
                   const oops::GeometryData &,
                   const oops::Variables &,
                   const Parameters_ &,
                   const State_ &,
                   const State_ &,
                   const double & adjointTolerance = -1.0);
  ~SaberOuterTBlock() {}

  const oops::GeometryData & innerGeometryData() {return saberOuterBlock_->innerGeometryData();}
  const oops::Variables & innerVars() {return saberOuterBlock_->innerVars();}

  void multiply(atlas::FieldSet &) const;
  void multiplyAD(atlas::FieldSet &) const;
  void calibrationInverseMultiply(atlas::FieldSet &) const;

 private:
  void print(std::ostream &);
  std::unique_ptr<SaberOuterBlockBase> saberOuterBlock_;
};

// ----------------------------------------------------------------------------

template <typename MODEL>
SaberOuterTBlock<MODEL>::SaberOuterTBlock(const Geometry_ & geom,
                                          const oops::GeometryData & geometryData,
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

  // Create outer block
  saberOuterBlock_.reset(SaberOuterBlockFactory::create(
                         geometryData,
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

  if (adjointTolerance >= 0.0) {
    // Adjoint test

    // Variables sizes
    std::vector<size_t> innerVariableSizes = geom.variableSizes(innerVars);

    // Create random inner FieldSet
    atlas::FieldSet innerFset = util::createRandomFieldSet(geometryData.functionSpace(),
                                                           innerVariableSizes,
                                                           innerVars);

    // Copy inner FieldSet
    atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

    // Variables sizes
    std::vector<size_t> outerVariableSizes = geom.variableSizes(outerVars);

    // Create random outer FieldSet
    atlas::FieldSet outerFset =
      util::createRandomFieldSet(geometryData.functionSpace(),
                                 outerVariableSizes,
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
