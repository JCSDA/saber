/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_HYDROBALSABERBLOCK_H_
#define SABER_VADER_HYDROBALSABERBLOCK_H_

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

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
class HydroBalSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydroBalSaberBlockParameters, SaberBlockParametersBase)
 public:
};

// -----------------------------------------------------------------------------
// This saber block is here
//
// -----------------------------------------------------------------------------
template <typename MODEL>
class HydroBalSaberBlock : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>             Geometry_;
  typedef oops::Increment<MODEL>            Increment_;
  typedef oops::State<MODEL>                State_;

 public:
  static const std::string classname() {return "saber::HydroBalSaberBlock";}

  typedef HydroBalSaberBlockParameters<MODEL> Parameters_;

  HydroBalSaberBlock(const Geometry_ &,
         const Parameters_ &,
         const State_ &,
         const State_ &);
  virtual ~HydroBalSaberBlock();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  oops::Variables inputVars_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
HydroBalSaberBlock<MODEL>::HydroBalSaberBlock(const Geometry_ & resol,
                      const Parameters_ & params,
                      const State_ & xb,
                      const State_ & fg)
  : SaberBlockBase<MODEL>(params),
    inputVars_(params.inputVars.value()),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydroBalSaberBlock starting" << std::endl;

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
    "virtual_potential_temperature"
  };

  std::vector<std::string> requiredGeometryVariables{"height_levels"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  // Use meta data to see if they are populated with actual data.
  for (auto & s : requiredStateVariables) {
    if (!xb.variables().has(s)) {
      oops::Log::info() << "HydroBalSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s]);
  }

  // check how virtual potential temperature is calculated.
  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(resol.extraFields()[s]);
  }

  for (auto & fld : augmentedStateFieldSet_) {
    double zz(0.0);
    auto view1 = atlas::array::make_view<double, 2>(fld);
    for (atlas::idx_t jnode = 0; jnode < fld.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < fld.shape(1); ++jlevel) {
        zz += view1(jnode, jlevel) * view1(jnode, jlevel);
      }
    }
    std::cout << "norm state fld :: " << fld.name() << " " << zz << std::endl;
  }

  oops::Log::trace() << classname() << "::HydroBalSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
HydroBalSaberBlock<MODEL>::~HydroBalSaberBlock() {
  oops::Log::trace() << classname() << "::~HydroBalSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~HydroBalSaberBlock");
  oops::Log::trace() << classname() << "::~HydroBalSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("HydroBalSaberBlock<MODEL>::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  for (auto & fld : fset) {
    double zz(0.0);
    auto view1 = atlas::array::make_view<double, 2>(fld);
    for (atlas::idx_t jnode = 0; jnode < fld.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < fld.shape(1); ++jlevel) {
        zz += view1(jnode, jlevel)*view1(jnode, jlevel);
      }
    }
    std::cout << "norm state inc before fld :: " << fld.name() << " " << zz << std::endl;
  }


  mo::hexner2ThetavTL(fset, augmentedStateFieldSet_);

  for (auto & fld : fset) {
    double zz(0.0);
    auto view1 = atlas::array::make_view<double, 2>(fld);
    for (atlas::idx_t jnode = 0; jnode < fld.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < fld.shape(1); ++jlevel) {
        zz += view1(jnode, jlevel)*view1(jnode, jlevel);
      }
    }
    std::cout << "norm state inc after fld :: " << fld.name() << " " << zz << std::endl;
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  mo::thetavP2HexnerTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::hexner2ThetavAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  mo::thetavP2HexnerAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HydroBalSaberBlock<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_HYDROBALSABERBLOCK_H_
