/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydroBalSaberBlock.h"

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

#include "oops/base/Variables.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<HydroBalSaberBlock>
       makerHydroBalSaberBlock_("mo_hydro_bal");

// -----------------------------------------------------------------------------

HydroBalSaberBlock::HydroBalSaberBlock(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & inputFunctionSpace,
               const atlas::FieldSet & inputExtraFields,
               const std::vector<size_t> & inputVariableSizes,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & outputVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::HydroBalSaberBlock starting" << std::endl;

  // Deserialize configuration
  HydroBalSaberBlockParameters params;
  params.deserialize(conf);

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
    if (!xb.has_field(s)) {
      oops::Log::info() << "HydroBalSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  // check how virtual potential temperature is calculated.
  mo::evalAirTemperature(augmentedStateFieldSet_);
  mo::evalTotalMassMoistAir(augmentedStateFieldSet_);
  mo::evalSatVaporPressure(augmentedStateFieldSet_);
  mo::evalSatSpecificHumidity(augmentedStateFieldSet_);
  mo::evalSpecificHumidity(augmentedStateFieldSet_);
  mo::evalVirtualPotentialTemperature(augmentedStateFieldSet_);

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(inputExtraFields[s]);
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

HydroBalSaberBlock::~HydroBalSaberBlock() {
  oops::Log::trace() << classname() << "::~HydroBalSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~HydroBalSaberBlock");
  oops::Log::trace() << classname() << "::~HydroBalSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalSaberBlock::multiply(atlas::FieldSet & fset) const {
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

void HydroBalSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::hexner2ThetavAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalSaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  mo::thetavP2HexnerAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydroBalSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
