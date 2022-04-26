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

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/util/atlasFieldSet.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(StdDevParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> fileConfig{"file", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class StdDev : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::StdDev";}

  typedef StdDevParameters Parameters_;

  StdDev(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~StdDev();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<atlas::FieldSet> stdDevFieldSet_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
StdDev<MODEL>::StdDev(const Geometry_ & resol,
                      const StdDevParameters & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params)
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

  // Setup increment
  Increment_ stdDev(resol, activeVars, xb.validTime());
  stdDev.read(params.fileConfig.value());
  oops::Log::test() << "Norm of stddev: " << std::scientific
                    << std::setprecision(3) << stdDev.norm() << std::endl;

  // Increment_ to ATLAS fieldset
  stdDevFieldSet_.reset(new atlas::FieldSet());
  stdDev.setAtlas(stdDevFieldSet_.get());
  stdDev.toAtlas(stdDevFieldSet_.get());

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StdDev<MODEL>::~StdDev() {
  oops::Log::trace() << classname() << "::~StdDev starting" << std::endl;
  util::Timer timer(classname(), "~StdDev");
  oops::Log::trace() << classname() << "::~StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::randomize(atlas::FieldSet * incFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  this->multiply(incFieldSet);
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::multiply(atlas::FieldSet * incFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  multiplyAtlasFieldSet(incFieldSet, stdDevFieldSet_.get());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::inverseMultiply(atlas::FieldSet * incFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  divideAtlasFieldSet(incFieldSet, stdDevFieldSet_.get());
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::multiplyAD(atlas::FieldSet * incFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  this->multiply(incFieldSet);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::inverseMultiplyAD(atlas::FieldSet * incFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  this->inverseMultiply(incFieldSet);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StdDev<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_STDDEV_H_
