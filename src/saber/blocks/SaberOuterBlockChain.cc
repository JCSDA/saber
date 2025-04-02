/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>
#include <tuple>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberOuterBlockChain::SaberOuterBlockChain(
                     const oops::GeometryData & outerGeometryData,
                     const oops::Variables & outerVars,
                     oops::FieldSet4D & fset4dXb,
                     oops::FieldSet4D & fset4dFg,
                     const eckit::LocalConfiguration & covarConf,
                     const std::vector<SaberOuterBlockParametersWrapper> & params) {
  oops::Log::trace() << "SaberOuterBlockChain generic ctor starting" << std::endl;
  oops::Log::info() << "Info     : Creating outer blocks" << std::endl;

  // Note model data information is not passed to vader.
  // This could cause issues for blocks which depend on it.

  // TODO(AS): check whether covarConf needs to be passed to the blocks (ideally not)
  const eckit::LocalConfiguration outerBlockConf{covarConf};

  // Loop in reverse order
  for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
    boost::adaptors::reverse(params)) {
    // Initialize current outer geometry data
    const oops::GeometryData & currentOuterGeometryData = outerBlocks_.size() == 0 ?
                                     outerGeometryData : outerBlocks_.back()->innerGeometryData();

    // Initialize outer block
    const auto[saberOuterBlockParams,
               currentOuterVars,
               activeVars]
            = initBlock(saberOuterBlockParamWrapper,
                          outerBlockConf,
                          currentOuterGeometryData,
                          outerVars,
                          fset4dXb,
                          fset4dFg);

    // Check block doesn't expect model fields to be read as this is a generic ctor
    if (outerBlocks_.back()->getReadConfs().size() != 0) {
      throw eckit::UserError("The generic constructor of the SABER outer block chain "
                             "does not allow to read MODEL fields.", Here());
    }

    // Check block doesn't expect calibration, as this could be done with the standard ctor
    if (saberOuterBlockParams.doCalibration()) {
      throw eckit::UserError("The generic constructor of the SABER outer block chain "
                             "does not allow covariance calibration.", Here());
    }

    if (saberOuterBlockParams.doRead()) {
      // Read data
      oops::Log::info() << "Info     : Read data" << std::endl;
      outerBlocks_.back()->read();

      if (saberOuterBlockParams.forceWrite.value()) {
        // Write data
        oops::Log::info() << "Info     : Write data" << std::endl;
        outerBlocks_.back()->write();
      }
    }


    // Inner geometry data and variables & consistency check with active variables
    auto[innerGeometryData, innerVars] = getInnerObjects(activeVars, currentOuterVars);

    // Left inverse multiplication on xb and fg if inner and outer Geometry are different
    interpolateStates(saberOuterBlockParams,
                      currentOuterGeometryData,
                      innerGeometryData,
                      fset4dXb,
                      fset4dFg);

    // Adjoint and inverse tests
    testLastOuterBlock(covarConf,
                       saberOuterBlockParams,
                       currentOuterGeometryData,
                       currentOuterVars,
                       innerGeometryData,
                       innerVars,
                       activeVars);
  }
  oops::Log::trace() << "SaberOuterBlockChain generic ctor done" << std::endl;
}

// -----------------------------------------------------------------------------
std::tuple<const SaberBlockParametersBase&, oops::Variables, oops::Variables>
    SaberOuterBlockChain::initBlock(
            const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper,
            const eckit::LocalConfiguration & outerBlockConf,
            const oops::GeometryData & outerGeometryData,
            const oops::Variables & outerVars,
            oops::FieldSet4D & fset4dXb,
            oops::FieldSet4D & fset4dFg) {
  // Initialize current outer variables and outer geometry data
  const oops::Variables & currentOuterVars = outerBlocks_.size() == 0 ?
                                   outerVars : outerBlocks_.back()->innerVars();

  // Get outer block parameters
  const SaberBlockParametersBase & saberOuterBlockParams =
    saberOuterBlockParamWrapper.saberOuterBlockParameters;
  oops::Log::info() << "Info     : Creating outer block: "
                    << saberOuterBlockParams.saberBlockName.value() << std::endl;

  // Get active variables
  const oops::Variables activeVars = getActiveVars(saberOuterBlockParams, currentOuterVars);

  // Get required variables in xb, fg if needed
  const oops::Variables mandatoryStateVars = saberOuterBlockParams.mandatoryStateVars();
  if (!(mandatoryStateVars <= fset4dXb.variables())) {
    oops::Log::info() << "Info     : Calling vader to populate trajectory variables for the "
                      << saberOuterBlockParams.saberBlockName.value() << std::endl;
    eckit::LocalConfiguration vaderCookbookConfig, vaderConfig;
    for (const auto & entry : saberDefaultCookbook) {
      vaderCookbookConfig.set(entry.first.name(), entry.second);
    }
    vaderConfig.set(vader::configCookbookKey, vaderCookbookConfig);
    vader::VaderParameters vaderParams;
    vader::Vader vader(vaderParams, vaderConfig);
    oops::Variables varsToPopulateXb(mandatoryStateVars);
    oops::Variables varsPopulatedXb = vader.changeVar(fset4dXb[0].fieldSet(),
                                                      varsToPopulateXb);
    if (varsToPopulateXb.size() != 0) {
      std::stringstream errorMsg;
      errorMsg << "Vader could not produce the requested variables "
               << varsToPopulateXb.variables()
               << " in block " << saberOuterBlockParams.saberBlockName.value()
               << std::endl;
      throw eckit::Exception(errorMsg.str(), Here());
    }
    oops::Variables varsToPopulateFg(mandatoryStateVars);
    oops::Variables varsPopulatedFg = vader.changeVar(fset4dFg[0].fieldSet(),
                                                      varsToPopulateFg);
    if (varsToPopulateFg.size() != 0) {
      std::stringstream errorMsg;
      errorMsg << "Vader could not produce the requested variables "
               << varsToPopulateFg.variables()
               << " in block " << saberOuterBlockParams.saberBlockName.value()
               << std::endl;
      throw eckit::Exception(errorMsg.str(), Here());
    }
  }

  // Create outer block
  outerBlocks_.emplace_back(SaberOuterBlockFactory::create(
                                               outerGeometryData,
                                               currentOuterVars,
                                               outerBlockConf,
                                               saberOuterBlockParams,
                                               fset4dXb[0],
                                               fset4dFg[0]));

  return std::tuple<const SaberBlockParametersBase&, oops::Variables, oops::Variables>(
              saberOuterBlockParams, currentOuterVars, activeVars);
}

// -----------------------------------------------------------------------------

std::tuple<const oops::GeometryData &, const oops::Variables &>
    SaberOuterBlockChain::getInnerObjects(const oops::Variables & activeVars,
                                          const oops::Variables & outerVars) const {
  // Inner variables and inner geometry data
  const oops::Variables & innerVars = outerBlocks_.back()->innerVars();
  const oops::GeometryData & innerGeometryData = outerBlocks_.back()->innerGeometryData();

  // Check that active variables are present in either inner or outer variables, or both
  for (const auto & var : activeVars) {
    if (!(innerVars.has(var) || outerVars.has(var))) {
      throw eckit::UserError("Active variable " + var.name() + " is not present in inner "
                             "or outer variables", Here());
    }
  }

  return std::tuple<const oops::GeometryData &, const oops::Variables &>(
              innerGeometryData, innerVars);
}

// -----------------------------------------------------------------------------

void SaberOuterBlockChain::interpolateStates(
        const SaberBlockParametersBase & saberOuterBlockParams,
        const oops::GeometryData & outerGeometryData,
        const oops::GeometryData & innerGeometryData,
        oops::FieldSet4D & fset4dXb,
        oops::FieldSet4D & fset4dFg) const {
  // Left inverse multiplication on xb and fg if inner and outer Geometry are different
  if (util::getGridUid(innerGeometryData.functionSpace())
    != util::getGridUid(outerGeometryData.functionSpace())
    && saberOuterBlockParams.inverseVars.value().size() > 0) {
    oops::Log::info() << "Info     : Left inverse multiplication on xb and fg" << std::endl;

    // Apply left inverse
    for (size_t itime = 0; itime < fset4dXb.size(); ++itime) {
      outerBlocks_.back()->leftInverseMultiply(fset4dXb[itime]);
      outerBlocks_.back()->leftInverseMultiply(fset4dFg[itime]);
    }
  }
}

// -----------------------------------------------------------------------------

void SaberOuterBlockChain::testLastOuterBlock(
                        const eckit::LocalConfiguration & covarConf,
                        const SaberBlockParametersBase & saberOuterBlockParams,
                        const oops::GeometryData & outerGeometryData,
                        const oops::Variables & outerVars,
                        const oops::GeometryData & innerGeometryData,
                        const oops::Variables & innerVars,
                        const oops::Variables & activeVars) const {
  // Get intersection of active variables and outer/inner variables
  oops::Variables activeOuterVars = outerVars;
  activeOuterVars.intersection(activeVars);
  oops::Variables activeInnerVars = innerVars;
  activeInnerVars.intersection(activeVars);

  // Adjoint test
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberOuterBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    outerBlocks_.back()->adjointTest(outerGeometryData,
                                     activeOuterVars,
                                     innerGeometryData,
                                     activeInnerVars,
                                     localAdjointTolerance);
  }

  // Inverse test
  const bool skipInverseTest = saberOuterBlockParams.skipInverseTest.value();
  if (covarConf.getBool("inverse test", false)) {
    oops::Log::info() << "Info     : Inverse test" << std::endl;
    if (skipInverseTest) {
      oops::Log::test() << "skipping inverse test for block "
                        << outerBlocks_.back()->blockName() << std::endl;
    } else {
      // Get inner and outer tolerances
      const double innerInverseTolerance = saberOuterBlockParams.innerInverseTolerance.value()
        .get_value_or(covarConf.getDouble("inverse tolerance"));
      const double outerInverseTolerance = saberOuterBlockParams.outerInverseTolerance.value()
        .get_value_or(covarConf.getDouble("inverse tolerance"));

      // Get inner and outer variables to compare
      oops::Variables innerVarsToCompare = saberOuterBlockParams.innerVariables.value()
        .get_value_or(activeInnerVars);
      oops::Variables outerVarsToCompare = saberOuterBlockParams.outerVariables.value()
        .get_value_or(activeOuterVars);

      // Run test
      outerBlocks_.back()->inverseTest(innerGeometryData,
                                       activeInnerVars,
                                       outerGeometryData,
                                       activeOuterVars,
                                       innerVarsToCompare,
                                       outerVarsToCompare,
                                       innerInverseTolerance,
                                       outerInverseTolerance);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
