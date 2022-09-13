/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_
#define SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "mo/common_varchange.h"
#include "mo/control2analysis_linearvarchange.h"
#include "mo/control2analysis_varchange.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/HydrostaticExnerParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
class HydrostaticExnerSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerSaberBlockParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<hydrostaticexnerParameters<MODEL>>
    hydrostaticexnerParams{"covariance data", this};
};

// -----------------------------------------------------------------------------
// This saber block is here to do 3 jobs
// 1) the vertical regression on geostrophic pressure
// 2) summing the result with unbalanced pressure to create hydrostatic_pressure
// 3) converting hydrostatic pressure to exner pressure.
// -----------------------------------------------------------------------------
template <typename MODEL>
class HydrostaticExnerSaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::HydrostaticExnerSaberBlock";}

  typedef HydrostaticExnerSaberBlockParameters<MODEL> Parameters_;

  HydrostaticExnerSaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~HydrostaticExnerSaberBlock();

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
HydrostaticExnerSaberBlock<MODEL>::HydrostaticExnerSaberBlock(const Geometry_ & resol,
                      const Parameters_ & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params),
    inputVars_(params.inputVars.value()),
    covFieldSet_(createGpRegressionStats(resol,
                                         inputVars_, params.hydrostaticexnerParams.value())),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydrostaticExnerSaberBlock starting" << std::endl;

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
    "air_pressure_levels_minus_one",
    "exner_levels_minus_one",
    "exner",
    "potential_temperature",
    "air_pressure_levels",
    "air_pressure",
    "m_v", "m_ci", "m_cl", "m_r",  // mixing ratios from file
    "m_t",  //  to be populated in evalTotalMassMoistAir
    "svp", "dlsvpdT",  //  to be populated in evalSatVaporPressure
    "qsat",  // to be populated in evalSatSpecificHumidity
    "specific_humidity",  //  to be populated in evalSpecificHumidity
    "virtual_potential_temperature",
    "hydrostatic_exner_levels", "hydrostatic_pressure_levels"
     };

  std::vector<std::string> requiredGeometryVariables{"height_levels"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::info() << "HydrostaticExnerSaberBlock variable " << s <<
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

  // we will need geometry here for height variables.
  mo::evalAirPressureLevels(augmentedStateFieldSet_);
  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);
  mo::evalHydrostaticExnerLevels(augmentedStateFieldSet_);
  mo::evalHydrostaticPressureLevels(augmentedStateFieldSet_);


  // Need to setup derived state fields that we need.
  std::vector<std::string> requiredCovarianceVariables{
    "vertical_regression_matrices", "interpolation_weights"};

  for (const auto & s : requiredCovarianceVariables) {
    augmentedStateFieldSet_.add(covFieldSet_[s]);
  }

  oops::Log::trace() << classname() << "::HydrostaticExnerSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
HydrostaticExnerSaberBlock<MODEL>::~HydrostaticExnerSaberBlock() {
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~HydrostaticExnerSaberBlock");
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("HydrostaticExnerSaberBlock<MODEL>::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalHydrostaticPressureTL(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticExnerTL(fset, augmentedStateFieldSet_);
  const auto hydrostaticPressureView =
      atlas::array::make_view<const double, 2>(fset["hydrostatic_pressure_levels"]);
  auto airPressureView =
      atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  airPressureView.assign(hydrostaticPressureView);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  auto airPressureView =
      atlas::array::make_view<double, 2>(fset["air_pressure_levels"]);
  auto hydrostaticPressureView =
      atlas::array::make_view<double, 2>(fset["hydrostatic_pressure_levels"]);

  for (atlas::idx_t jn = 0; jn < fset["hydrostatic_exner_levels"].shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < fset["hydrostatic_exner_levels"].shape(1); ++jl) {
      hydrostaticPressureView(jn, jl) += airPressureView(jn, jl);
      airPressureView(jn, jl) = 0.0;
    }
  }
  mo::evalHydrostaticExnerAD(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticPressureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiplyAD not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydrostaticExnerSaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_
