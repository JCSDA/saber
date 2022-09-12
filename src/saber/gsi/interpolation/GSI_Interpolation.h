/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
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

class InterpolationParameters : public InterpolationImplParameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, InterpolationImplParameters)
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class Interpolation : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;

 public:
  static const std::string classname() {return "saber::gsi::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const Geometry_ &, const Parameters_ &, const State_ &, const State_ &);
  virtual ~Interpolation();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  // GSI Interpolation implementation
  InterpolationImpl interpolationImpl_;
};

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
Interpolation<MODEL>::Interpolation(const Geometry_ & geom, const Parameters_ & params,
                                    const State_ & xb, const State_ & fg)
  : SaberBlockBase<MODEL>(params),
interpolationImpl_(geom.getComm(), geom.functionSpace(), params,
  params.inputVars.value().variables())
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;
  util::Timer timer(classname(), "Interpolation");
  // Assert that there is no variable change in this block
  ASSERT(params.inputVars.value() == params.outputVars.value());
  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
Interpolation<MODEL>::~Interpolation() {
  oops::Log::trace() << classname() << "::~Interpolation starting" << std::endl;
  util::Timer timer(classname(), "~Interpolation");
  oops::Log::trace() << classname() << "::~Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::randomize(atlas::FieldSet & fset) const {
  ABORT(classname() + "randomize: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  interpolationImpl_.multiply(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiply: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  interpolationImpl_.multiplyAD(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiplyAD: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
