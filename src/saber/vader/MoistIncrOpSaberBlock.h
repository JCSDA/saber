/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_MOISTINCROPSABERBLOCK_H_
#define SABER_VADER_MOISTINCROPSABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/functions.h"
#include "mo/model2geovals_varchange.h"

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

class MoistIncrOpSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistIncrOpSaberBlockParameters, SaberBlockParametersBase)
 public:
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class MoistIncrOpSaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::AirTemperatureSaberBlock";}

  typedef MoistIncrOpSaberBlockParameters Parameters_;

  MoistIncrOpSaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~MoistIncrOpSaberBlock();

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
MoistIncrOpSaberBlock<MODEL>::MoistIncrOpSaberBlock(const Geometry_ &,  // resol,
                      const MoistIncrOpSaberBlockParameters & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistIncrOpSaberBlock starting" << std::endl;

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
  std::vector<std::string> requiredStateVariables{
    "potential_temperature",  // from file
    "exner",  // on theta levels from file ("exner_levels_minus_one" is on rho levels)
    "air_pressure",  // on theta levels from file ("air_pressure_levels_minus_one" is on rho levels)
    "air_temperature",  // to be populated in evalAirTemperature
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  // to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  // to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  // to be populated in evalSpecificHumidity
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in evalMassCloudLiquid
    "mass_content_of_cloud_ice_in_atmosphere_layer",  // to be populated in evalMassCloudIce
    "qrain",  // to be populated in evalMassRain
    "rht",  // to be populated in evalTotalRelativeHumidity
    "liquid_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "ice_cloud_volume_fraction_in_atmosphere_layer",  // from file
    "cleff", "cfeff"  // to be populated in getMIOFields
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::info() << "MoistIncrOpSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalMassCloudLiquid(augmentedStateFieldSet_);
  mo::evalMassCloudIce(augmentedStateFieldSet_);
  mo::evalMassRain(augmentedStateFieldSet_);
  mo::evalTotalRelativeHumidity(augmentedStateFieldSet_);
  mo::functions::getMIOFields(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::MoistIncrOpSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
MoistIncrOpSaberBlock<MODEL>::~MoistIncrOpSaberBlock() {
  oops::Log::trace() << classname() << "::~MoistIncrOpSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~MoistIncrOpSaberBlock");
  oops::Log::trace() << classname() << "::~MoistIncrOpSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("MoistIncrOpSaberBlock<MODEL>::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::qtTemperature2qqclqcfTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  mo::qqclqcf2qtTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::qtTemperature2qqclqcfAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  mo::qqclqcf2qtAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistIncrOpSaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_MOISTINCROPSABERBLOCK_H_
