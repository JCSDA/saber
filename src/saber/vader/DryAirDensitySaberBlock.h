/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_DRYAIRDENSITYSABERBLOCK_H_
#define SABER_VADER_DRYAIRDENSITYSABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
// -----------------------------------------------------------------------------

class DryAirDensitySaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensitySaberBlockParameters, SaberBlockParametersBase)
 public:
};

template <typename MODEL>
class DryAirDensitySaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::DryAirDensitySaberBlock";}

  typedef DryAirDensitySaberBlockParameters Parameters_;

  DryAirDensitySaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);

  virtual ~DryAirDensitySaberBlock();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  atlas::FieldSet augmentedStateFieldSet_;
};  // class definition DryAirDensitySaberBlock

// -----------------------------------------------------------------------------

template<typename MODEL>
DryAirDensitySaberBlock<MODEL>::DryAirDensitySaberBlock(const Geometry_ & resol,
                      const DryAirDensitySaberBlockParameters & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::DryAirDensitySaberBlock starting" << std::endl;

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
  std::vector<std::string> requiredStateVariables{ "exner_levels_minus_one",
                                                   "potential_temperature",
                                                   "exner",
                                                   "air_pressure_levels_minus_one",
                                                   "air_temperature",
                                                   "dry_air_density_levels_minus_one"};

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::error() << "::DryAirDensitySaberBlock variable " << s <<
                            "is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(resol.extraFields()[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalDryAirDensity(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::DryAirDensitySaberBlock done" << std::endl;
}  // class declaration DryAirDensitySaberBlock

template<typename MODEL>
DryAirDensitySaberBlock<MODEL>::~DryAirDensitySaberBlock() {
  oops::Log::trace() << classname() << "::~DryAirDensitySaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~DryAirDensitySaberBlock");
  oops::Log::trace() << classname() << "::~DryAirDensitySaberBlock done" << std::endl;
}  // ~DryAirDensitySaberBlock

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("DryAirDensitySaberBlock<MODEL>::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}  // randomize

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalDryAirDensityTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}  // multiply

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}  // inverseMultiply

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalDryAirDensityAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}  // multiplyAD

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiplyAD not meaningful so fieldset unchanged"
                    << std::endl;
}  // inverseMultiplyAD

template<typename MODEL>
void DryAirDensitySaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

}  // namespace saber

#endif  // SABER_VADER_DRYAIRDENSITYSABERBLOCK_H_
