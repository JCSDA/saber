/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_AIRTEMPERATURESABERBLOCK_H_
#define SABER_VADER_AIRTEMPERATURESABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "mo/control2analysis_linearvarchange.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class AirTemperatureSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(AirTemperatureSaberBlockParameters, SaberBlockParametersBase)
 public:
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class AirTemperatureSaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::AirTemperatureSaberBlock";}

  typedef AirTemperatureSaberBlockParameters Parameters_;

  AirTemperatureSaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~AirTemperatureSaberBlock();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
AirTemperatureSaberBlock<MODEL>::AirTemperatureSaberBlock(const Geometry_ & resol,
                      const AirTemperatureSaberBlockParameters & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::AirTemperatureSaberBlock starting" << std::endl;

  // Setup and check input/output variables
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

  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredStateVariables{"exner_levels_minus_one",
                                                  "potential_temperature"};

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::info() << "AirTemperatureSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(resol.extraFields()[s]);
  }

  oops::Log::trace() << classname() << "::AirTemperatureSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
AirTemperatureSaberBlock<MODEL>::~AirTemperatureSaberBlock() {
  oops::Log::trace() << classname() << "::~AirTemperatureSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~AirTemperatureSaberBlock");
  oops::Log::trace() << classname() << "::~AirTemperatureSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  ABORT("AirTemperatureSaberBlock<MODEL>::randomize: not implemented");
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalAirTemperatureTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;

  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalAirTemperatureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  this->inverseMultiply(fset);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void AirTemperatureSaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_AIRTEMPERATURESABERBLOCK_H_
