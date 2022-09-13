/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/HydrostaticExnerSaberBlock.h"

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

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/HydrostaticExnerParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberBlockMaker<HydrostaticExnerSaberBlock>
       makerHydrostaticExnerSaberBlock_("mo_hydrostatic_exner");

// -----------------------------------------------------------------------------

HydrostaticExnerSaberBlock::HydrostaticExnerSaberBlock(const eckit::mpi::Comm & comm,
                                                       const atlas::FunctionSpace & functionSpace,
                                                       const atlas::FieldSet & extraFields,
                                                       const std::vector<size_t> & variableSizes,
                                                       const Parameters_ & params,
                                                       const atlas::FieldSet & xb,
                                                       const atlas::FieldSet & fg,
                                                       const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params),
    inputVars_(params.inputVars.value()),
    covFieldSet_(createGpRegressionStats(functionSpace, extraFields, variableSizes,
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
    if (!xb.has_field(s)) {
      oops::Log::info() << "HydrostaticExnerSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    augmentedStateFieldSet_.add(extraFields[s]);
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

HydrostaticExnerSaberBlock::~HydrostaticExnerSaberBlock() {
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~HydrostaticExnerSaberBlock");
  oops::Log::trace() << classname() << "::~HydrostaticExnerSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  throw eckit::NotImplemented("HydrostaticExnerSaberBlock::randomize", Here());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalHydrostaticPressureTL(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticExnerTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalHydrostaticExnerAD(fset, augmentedStateFieldSet_);
  mo::evalHydrostaticPressureAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::inverseMultiplyAD not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void HydrostaticExnerSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
