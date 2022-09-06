/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_MOISTURECONTROLSABERBLOCK_H_
#define SABER_VADER_MOISTURECONTROLSABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/MoistureControlParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
class MoistureControlSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistureControlSaberBlockParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<moisturecontrolParameters<MODEL>>
    moisturecontrolParams{"covariance data", this};
};

// -----------------------------------------------------------------------------
// This saber block is here
//
// -----------------------------------------------------------------------------
template <typename MODEL>
class MoistureControlSaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::MoistureControlSaberBlock";}

  typedef MoistureControlSaberBlockParameters<MODEL> Parameters_;

  MoistureControlSaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~MoistureControlSaberBlock();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  oops::Variables inputVars_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
MoistureControlSaberBlock<MODEL>::MoistureControlSaberBlock(const Geometry_ & resol,
                      const Parameters_ & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params),
    inputVars_(params.inputVars.value()),
    covFieldSet_(createMuStats(resol, params.moisturecontrolParams.value())),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistureControlSaberBlock starting" << std::endl;

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

  std::vector<std::string> requiredStateVariables{
    "air_temperature",
    "air_pressure",
    "potential_temperature",   // from file
    "exner",  // from file on theta levels ("exner_levels_minus_one" is on rho levels)
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  //  to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  //  to be populated in evalSpecificHumidity
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
      // to be populated in evalMassCloudLiquid
    "mass_content_of_cloud_ice_in_atmosphere_layer",  // to be populated in evalMassCloudIce
    "qrain",  // to be populated in evalMassRain
    "qt",  // to be populated in qqclqcf2qtTL
    "rht",  // to be populated in evalTotalRelativeHumidity
    "muA", "muH1",  // to be populated in function call from CovarianceStatisticsUtils.h
    "muRow1Column1", "muRow1Column2",  // to be populated in evalMoistureControlDependencies
    "muRow2Column1", "muRow2Column2",  //   ""
    "muRecipDeterminant"  //   ""
  };

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::info() << "MoistureControlSaberBlock variable " << s <<
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
  mo::qqclqcf2qt(augmentedStateFieldSet_);
  mo::evalTotalRelativeHumidity(augmentedStateFieldSet_);

  // populate "muA" and "muH1"
  for (auto & covFld : covFieldSet_) {
    populateInterpMuStats(augmentedStateFieldSet_,
                          covFld);
  }

  // populate "specific moisture control dependencies"
  mo::evalMoistureControlDependencies(augmentedStateFieldSet_);

  oops::Log::trace() << classname() << "::MoistureControlSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
MoistureControlSaberBlock<MODEL>::~MoistureControlSaberBlock() {
  oops::Log::trace() << classname() << "::~MoistureControlSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~MoistureControlSaberBlock");
  oops::Log::trace() << classname() << "::~MoistureControlSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("MoistureControlSaberBlock<MODEL>::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalQtThetaTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  mo::evalMuThetavTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalQtThetaAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  mo::evalMuThetavAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MoistureControlSaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_MOISTURECONTROLSABERBLOCK_H_

