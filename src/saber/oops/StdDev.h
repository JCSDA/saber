/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_STDDEV_H_
#define SABER_OOPS_STDDEV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(StdDevParameters, SaberBlockParametersBase)
};

// -----------------------------------------------------------------------------

class StdDev : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::StdDev";}

  typedef StdDevParameters Parameters_;

  StdDev(const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const Parameters_ &,
         const atlas::FieldSet &,
         const atlas::FieldSet &);
  virtual ~StdDev();

  void initialize(const std::vector<atlas::FieldSet> &) override;

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  atlas::FieldSet stdDevFset_;
};

// -----------------------------------------------------------------------------

StdDev::StdDev(const atlas::FunctionSpace & functionSpace,
               const atlas::FieldSet & extraFields,
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg)
  : SaberBlockBase(params), stdDevFset_()
{
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

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

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

StdDev::~StdDev() {
  oops::Log::trace() << classname() << "::~StdDev starting" << std::endl;
  util::Timer timer(classname(), "~StdDev");
  oops::Log::trace() << classname() << "::~StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::initialize(const std::vector<atlas::FieldSet> & fsetVec)  {
  oops::Log::trace() << classname() << "::initialize starting" << std::endl;

  // Check fsetVec size
  std::cout << "fsetVec.size() " << fsetVec.size() << std::endl;
  ASSERT(fsetVec.size() == 1);

  // Check FieldSet name
  std::cout << "fsetVec[0].name() " << fsetVec[0].name() << std::endl;
  ASSERT(fsetVec[0].name() == "StdDev");

  // Copy std-dev fields
  stdDevFset_.clear();
  for (const auto & field : fsetVec[0]) {
      stdDevFset_.add(field);
  }

  oops::Log::trace() << classname() << "::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  ABORT("StdDev::randomize: not implemented");
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::FieldSetMultiply(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  util::FieldSetDivide(fset, stdDevFset_);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  this->inverseMultiply(fset);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_STDDEV_H_
