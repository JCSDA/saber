/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/MoistureControlSaberBlock.h"

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
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"
#include "saber/vader/MoistureControlParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<MoistureControlSaberBlock> 
       makerMoistureControlBlock_("mo_moisture_control");

// -----------------------------------------------------------------------------

MoistureControlSaberBlock::MoistureControlSaberBlock(const eckit::mpi::Comm & comm,
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
  : SaberOuterBlockBase(conf), augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::MoistureControlSaberBlock starting" << std::endl;

  // Deserialize configuration
  MoistureControlSaberBlockParameters params;
  params.deserialize(conf);

  // Covariance FieldSet
  covFieldSet_ = createMuStats(inputExtraFields, params.moisturecontrolParams.value());

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
    if (!xb.has_field(s)) {
      oops::Log::info() << "MoistureControlSaberBlock variable " << s <<
                           " is not part of state object." << std::endl;
    }
  }

  augmentedStateFieldSet_.clear();
  for (const auto & s : requiredStateVariables) {
    augmentedStateFieldSet_.add(xb[s]);
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

MoistureControlSaberBlock::~MoistureControlSaberBlock() {
  oops::Log::trace() << classname() << "::~MoistureControlSaberBlock starting" << std::endl;
  util::Timer timer(classname(), "~MoistureControlSaberBlock");
  oops::Log::trace() << classname() << "::~MoistureControlSaberBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControlSaberBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  mo::evalQtThetaTL(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControlSaberBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  mo::evalQtThetaAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControlSaberBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  mo::evalMuThetavAD(fset, augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void MoistureControlSaberBlock::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
