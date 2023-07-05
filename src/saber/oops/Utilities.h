/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/range/adaptor/reversed.hpp>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelData.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/SaberBlockChain.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {

oops::Variables getActiveVars(const SaberBlockParametersBase &,
                              const oops::Variables &);

// -----------------------------------------------------------------------------

template<typename MODEL>
eckit::LocalConfiguration readEnsemble(const oops::Geometry<MODEL> & geom,
                                       const oops::Variables & vars,
                                       const oops::State<MODEL> & xb,
                                       const oops::State<MODEL> & fg,
                                       const eckit::LocalConfiguration & inputConf,
                                       const bool & iterativeEnsembleLoading,
                                       std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Create output configuration
  eckit::LocalConfiguration outputConf;

  // Fill output configuration and set ensemble size
  size_t nens = 0;
  size_t ensembleFound = 0;

  // Ensemble of states, perturbation using the mean
  oops::IncrementEnsembleFromStatesParameters<MODEL> ensembleParams;
  eckit::LocalConfiguration ensembleConf;
  if (inputConf.has("ensemble")) {
    ensembleConf = inputConf.getSubConfiguration("ensemble");
    ensembleParams.deserialize(ensembleConf);
    nens = ensembleParams.states.size();
    outputConf.set("ensemble", ensembleConf);
    ++ensembleFound;
  }

  // Increment ensemble from increments on disk
  oops::IncrementEnsembleParameters<MODEL> ensemblePertParams;
  eckit::LocalConfiguration ensemblePert;
  if (inputConf.has("ensemble pert")) {
    ensemblePert = inputConf.getSubConfiguration("ensemble pert");
    ensemblePertParams.deserialize(ensemblePert);
    nens = ensemblePertParams.size();
    outputConf.set("ensemble pert", ensemblePert);
    ++ensembleFound;
  }

  // Increment ensemble from difference of two states
  oops::StateEnsembleParameters<MODEL> ensembleBaseParams;
  oops::StateEnsembleParameters<MODEL> ensemblePairsParams;
  eckit::LocalConfiguration ensembleBase;
  eckit::LocalConfiguration ensemblePairs;
  if (inputConf.has("ensemble base") && inputConf.has("ensemble pairs")) {
    ensembleBase = inputConf.getSubConfiguration("ensemble base");
    ensemblePairs = inputConf.getSubConfiguration("ensemble pairs");
    ensembleBaseParams.deserialize(ensembleBase);
    ensemblePairsParams.deserialize(ensemblePairs);
    nens = ensembleBaseParams.size();
    outputConf.set("ensemble base", ensembleBase);
    outputConf.set("ensemble pairs", ensemblePairs);
    ++ensembleFound;
  }

  // Set ensemble size
  outputConf.set("ensemble size", nens);

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  if (!iterativeEnsembleLoading) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble pointer
    std::unique_ptr<oops::IncrementEnsemble<MODEL>> ensemble;

    // Ensemble of states, perturbation using the mean
    if (!ensembleConf.empty()) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;
      ensemble.reset(new oops::IncrementEnsemble<MODEL>(ensembleParams, geom, vars,
        xb.validTime()));
    }

    // Increment ensemble from increments on disk
    if (!ensemblePert.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
      ensemble.reset(new oops::IncrementEnsemble<MODEL>(geom, vars, ensemblePertParams));
    }

    // Increment ensemble from difference of two states
    if (!ensembleBase.empty() && !ensemblePairs.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from difference of two states"
                        << std::endl;
      ensemble.reset(new oops::IncrementEnsemble<MODEL>(geom, vars, ensembleBaseParams,
        ensemblePairsParams));
    }

    // Transform Increment into FieldSet
    for (unsigned int ie = 0; ie < nens; ++ie) {
      atlas::FieldSet fset = util::shareFields((*ensemble)[ie].fieldSet());
      fset.name() = "ensemble member";
      fsetEns.push_back(fset);
    }
  }

  // Return ensemble configuration for iterative ensemble loading
  return outputConf;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::Variables & vars,
                      const util::DateTime & date,
                      const eckit::LocalConfiguration & conf,
                      atlas::FieldSet & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create Increment
  oops::Increment<MODEL> dx(geom, vars, date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset = util::shareFields(dx.fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::Variables & vars,
                        const util::DateTime & date,
                        const eckit::LocalConfiguration & conf,
                        const size_t & ie,
                        atlas::FieldSet & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  if (conf.has("ensemble")) {
    // Ensemble of states passed as increments
    oops::StateEnsembleParameters<MODEL> states;
    states.deserialize(conf.getSubConfiguration("ensemble"));

    // Read state
    oops::State<MODEL> xx(geom, states.getStateParameters(ie));

    // Copy FieldSet
    fset = util::copyFieldSet(xx.fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from increments on disk
    oops::IncrementEnsembleParameters<MODEL> ensemblePertParams;
    ensemblePertParams.deserialize(conf.getSubConfiguration("ensemble pert"));

    // Read Increment
    oops::Increment<MODEL> dx(geom, vars, date);
    dx.read(ensemblePertParams.getIncrementParameters(ie));

    // Get FieldSet
    fset = util::copyFieldSet(dx.fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble base") && conf.has("ensemble pairs")) {
    // Increment ensemble from difference of two states
    oops::StateEnsembleParameters<MODEL> ensembleBaseParams;
    ensembleBaseParams.deserialize(conf.getSubConfiguration("ensemble base"));
    oops::StateEnsembleParameters<MODEL> ensemblePairsParams;
    ensemblePairsParams.deserialize(conf.getSubConfiguration("ensemble pairs"));

    // Read states
    oops::State<MODEL> xxBase(geom, ensembleBaseParams.getStateParameters(ie));
    oops::State<MODEL> xxPairs(geom, ensemblePairsParams.getStateParameters(ie));

    // Compute difference
    oops::Increment<MODEL> dx(geom, vars, date);
    dx.diff(xxPairs, xxBase);

    // Get FieldSet
    fset = util::copyFieldSet(dx.fieldSet());

    ++ensembleFound;
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound <= 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void writeEnsemble(const oops::Geometry<MODEL> & geom,
                   const oops::Variables & vars,
                   const util::DateTime & date,
                   const OutputEnsembleParameters<MODEL> & outputEnsemble,
                   const eckit::LocalConfiguration & conf,
                   const bool & iterativeEnsembleLoading,
                   const SaberBlockChain & saberBlockChain,
                   const std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << "writeEnsemble starting" << std::endl;

  oops::Log::info() << "Info     : Write ensemble members" << std::endl;

  // Check whether geometry grid is similar to the last outer block inner geometry
  const bool useModelWriter = (util::getGridUid(geom.generic().functionSpace())
    == util::getGridUid(saberBlockChain.lastOuterBlock().innerGeometryData().functionSpace()));

  // Get ensemble size
  size_t ensembleSize = conf.getInt("ensemble size");

  // Estimate mean
  atlas::FieldSet mean;
  if (iterativeEnsembleLoading) {
    for (size_t ie = 0; ie < ensembleSize; ++ie) {
      // Read member
      atlas::FieldSet fset;
      readEnsembleMember(geom, vars, date, conf, ie, fset);

      // Update mean
      if (ie == 0) {
        mean = util::copyFieldSet(fset);
      } else {
        util::addFieldSets(mean, fset);
      }
    }

    // Normalize mean
    util::multiplyFieldSet(mean, 1.0/static_cast<double>(ensembleSize));
  }

  // Output parameters
  typename oops::Increment<MODEL>::WriteParameters_ writeParams
    = outputEnsemble.writeParams;
  const bool firstMemberOnly = outputEnsemble.firstMemberOnly.value();

  // Write first member only
  if (firstMemberOnly) {
    ensembleSize = 1;
  }

  for (size_t ie = 0; ie < ensembleSize; ++ie) {
    oops::Log::info() << "Info     : Write member " << ie << std::endl;

    // Increment pointer
    std::unique_ptr<oops::Increment<MODEL>> dx;

    // Get ensemble member
    if (iterativeEnsembleLoading) {
      // Read ensemble member
      dx.reset(new oops::Increment<MODEL>(geom, vars, date));
      readEnsembleMember(geom, vars, date, conf, ie, dx->fieldSet());

      // Remove mean
      util::subtractFieldSets(dx->fieldSet(), mean);

      // Apply outer blocks inverse
      saberBlockChain.leftInverseMultiply(dx->fieldSet());
    } else {
      // Copy member
      dx.reset(new oops::Increment<MODEL>(geom, vars, date));
      dx->fieldSet() = fsetEns[ie];
    }

    // ATLAS fieldset to Increment_
    dx->synchronizeFields();

    if (useModelWriter) {
      // Use model writer

      // Set member index
      writeParams.setMember(ie+1);

      // Write Increment
      dx->write(writeParams);
      oops::Log::test() << "Norm of ensemble member " << ie << ": " << dx->norm() << std::endl;
    } else {
      // Use generic ATLAS writer
      ABORT("generic output ensemble write not implemented yet");
    }
  }

  oops::Log::trace() << "writeEnsemble done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
std::unique_ptr<SaberBlockChain> ensembleBlockChain(const oops::Geometry<MODEL> & geom,
                                                    const oops::Variables & incVars,
                                                    const atlas::FieldSet & fsetXbIn,
                                                    const atlas::FieldSet & fsetFgIn,
                                                    const util::DateTime & validTimeOfXbFg,
                                                    const eckit::Configuration & conf) {
  oops::Log::trace() << "ensembleBlockChain starting" << std::endl;

  oops::Log::info() << "Info     : Create ensemble transform or localization blockchain"
                    << std::endl;

  // Extend backgroud and first guess with extra fields
  // TODO(Benjamin, Marek, Mayeul, ?)

  // Make copy of fieldsets that can be changed
  atlas::FieldSet fsetXb = util::copyFieldSet(fsetXbIn);
  atlas::FieldSet fsetFg = util::copyFieldSet(fsetFgIn);

  // Initialize geometry data
  std::vector<std::reference_wrapper<const oops::GeometryData>> outerGeometryData;
  outerGeometryData.push_back(geom.generic());

  // Intialize outer variables
  oops::Variables outerVars(incVars);

  // Initialize blockchain
  std::unique_ptr<SaberBlockChain> blockChain(new SaberBlockChain(incVars, fsetXb));

  // Initialize empty vector of FieldSets
  std::vector<atlas::FieldSet> emptyFsetEns;

  // Create covariance configuration
  eckit::LocalConfiguration covarConf;
  eckit::LocalConfiguration ensembleConf;
  ensembleConf.set("ensemble size", 0);
  covarConf.set("ensemble configuration", ensembleConf);
  covarConf.set("adjoint test", conf.getBool("adjoint test", false));
  covarConf.set("adjoint tolerance", conf.getDouble("adjoint tolerance", 1.0e-12));
  covarConf.set("inverse test", conf.getBool("inverse test", false));
  covarConf.set("inverse tolerance", conf.getDouble("inverse tolerance", 1.0e-12));
  covarConf.set("iterative ensemble loading", false);

  if (conf.has("saber outer blocks")) {
    // Build outer blocks successively
    std::vector<SaberOuterBlockParametersWrapper> saberOuterBlocksParams;
    for (const auto & saberOuterBlockConf : conf.getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper saberOuterBlockParamsWrapper;
      saberOuterBlockParamsWrapper.deserialize(saberOuterBlockConf);
      saberOuterBlocksParams.push_back(saberOuterBlockParamsWrapper);
      const SaberBlockParametersBase & saberOuterBlockParams =
        saberOuterBlockParamsWrapper.saberOuterBlockParameters;
      if (saberOuterBlockParams.doCalibration()) {
        ABORT("no calibration within ensemble transform / localization");
      }
    }
    buildOuterBlocks(geom,
                     outerGeometryData,
                     outerVars,
                     fsetXb,
                     fsetFg,
                     validTimeOfXbFg,
                     emptyFsetEns,
                     covarConf,
                     saberOuterBlocksParams,
                     *blockChain);
  }

  if (conf.has("saber central block")) {
    // Central block
    SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
    saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));
    const SaberBlockParametersBase & saberCentralBlockParams =
      saberCentralBlockParamsWrapper.saberCentralBlockParameters;
    if (saberCentralBlockParams.doCalibration()) {
      ABORT("no calibration within ensemble transform / localization");
    }
    buildCentralBlock(geom,
                      geom,
                      outerGeometryData.back().get(),
                      outerVars,
                      fsetXb,
                      fsetFg,
                      validTimeOfXbFg,
                      emptyFsetEns,
                      emptyFsetEns,
                      covarConf,
                      saberCentralBlockParamsWrapper,
                      *blockChain);
  }

  // Return blockchain
  return blockChain;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void buildOuterBlocks(const oops::Geometry<MODEL> & geom,
                      std::vector<std::reference_wrapper<const oops::GeometryData>> &
                        outerGeometryData,
                      oops::Variables & outerVars,
                      atlas::FieldSet & fsetXb,
                      atlas::FieldSet & fsetFg,
                      const util::DateTime & validTimeOfXbFg,
                      std::vector<atlas::FieldSet> & fsetEns,
                      const eckit::LocalConfiguration & covarConf,
                      const std::vector<SaberOuterBlockParametersWrapper> & saberOuterBlocksParams,
                      SaberBlockChain & saberBlockChain) {
  oops::Log::trace() << "buildOuterBlocks starting" << std::endl;

  oops::Log::info() << "Info     : Creating outer blocks" << std::endl;

  oops::ModelData<MODEL> modelData{geom};
  eckit::LocalConfiguration modelDataConf;
  modelDataConf.set("model data", modelData.modelData());  // Key here is vader::configModelVarsKey
  eckit::LocalConfiguration outerBlockConf{covarConf};
  outerBlockConf.set("vader", modelDataConf);

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Loop in reverse order
  for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
    boost::adaptors::reverse(saberOuterBlocksParams)) {
    // Get outer block parameters
    const SaberBlockParametersBase & saberOuterBlockParams =
      saberOuterBlockParamWrapper.saberOuterBlockParameters;
    oops::Log::info() << "Info     : Creating outer block: "
                      << saberOuterBlockParams.saberBlockName.value() << std::endl;

    // Get active variables
    oops::Variables activeVars = getActiveVars(saberOuterBlockParams, outerVars);

    // Create outer block
    saberBlockChain.outerBlocks().emplace_back(SaberOuterBlockFactory::create(
                                               outerGeometryData.back().get(),
                                               outerVars,
                                               outerBlockConf,
                                               saberOuterBlockParams,
                                               fsetXb,
                                               fsetFg,
                                               validTimeOfXbFg));

    // Read and add model fields
    saberBlockChain.lastOuterBlock().read(geom, outerVars, validTimeOfXbFg);

    if (saberOuterBlockParams.doCalibration()) {
      // Block calibration

      // Ensemble configuration
      eckit::LocalConfiguration ensembleConf
        = covarConf.getSubConfiguration("ensemble configuration");

      if (iterativeEnsembleLoading) {
        // Iterative calibration
         oops::Log::info() << "Info     : Iterative calibration" << std::endl;

        // Initialization
        saberBlockChain.lastOuterBlock().iterativeCalibrationInit();

        // Get ensemble size
        size_t nens = ensembleConf.getInt("ensemble size");

        for (size_t ie = 0; ie < nens; ++ie) {
          // Read ensemble member
          atlas::FieldSet fset;
          readEnsembleMember(geom,
                             saberBlockChain.incVars(),
                             validTimeOfXbFg,
                             ensembleConf,
                             ie,
                             fset);

          // Apply outer blocks inverse (except last)
          saberBlockChain.leftInverseMultiplyExceptLast(fset);

          // Use FieldSet in the central block
          oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
          saberBlockChain.lastOuterBlock().iterativeCalibrationUpdate(fset);
        }

        // Finalization
        oops::Log::info() << "Info     : Finalization" << std::endl;
        saberBlockChain.lastOuterBlock().iterativeCalibrationFinal();
      } else {
        // Direct calibration
        oops::Log::info() << "Info     : Direct calibration" << std::endl;
        saberBlockChain.lastOuterBlock().directCalibration(fsetEns);
      }

      // Write calibration data
      oops::Log::info() << "Info     : Write calibration data" << std::endl;
      saberBlockChain.lastOuterBlock().write(geom, outerVars, validTimeOfXbFg);
      saberBlockChain.lastOuterBlock().write();

      if (!iterativeEnsembleLoading) {
        // Left inverse multiplication on ensemble members
        oops::Log::info() << "Info     : Left inverse multiplication on ensemble members"
                          << std::endl;
        for (auto & fset : fsetEns) {
          if (saberBlockChain.lastOuterBlock().skipInverse()) {
            oops::Log::info()
                    << "Info     : Warning: left inverse multiplication skipped for block "
                    << saberBlockChain.lastOuterBlock().blockName() << std::endl;
          } else {
            saberBlockChain.lastOuterBlock().leftInverseMultiply(fset);
          }
        }
      }
    } else if (saberOuterBlockParams.doRead()) {
      // Read data
      oops::Log::info() << "Info     : Read data" << std::endl;
      saberBlockChain.lastOuterBlock().read();
    }

    // Inner geometry data and variables
    const oops::GeometryData & innerGeometryData =
      saberBlockChain.lastOuterBlock().innerGeometryData();
    const oops::Variables innerVars = saberBlockChain.lastOuterBlock().innerVars();

    // Check that active variables are present in either inner or outer variables, or both
    for (const auto & var : activeVars.variables()) {
      ASSERT(innerVars.has(var) || outerVars.has(var));
    }

    // Get intersection of active variables and outer/inner variables
    oops::Variables activeOuterVars = outerVars;
    activeOuterVars.intersection(activeVars);
    oops::Variables activeInnerVars = innerVars;
    activeInnerVars.intersection(activeVars);


    // Left inverse multiplication on xb and fg if inner and outer Geometry is different
    if (util::getGridUid(innerGeometryData.functionSpace())
      != util::getGridUid(outerGeometryData.back().get().functionSpace())
      && saberOuterBlockParams.inverseVars.value().size() > 0) {
      oops::Log::info() << "Info     : Left inverse multiplication on xb and fg" << std::endl;

      // Share fields pointers
      atlas::FieldSet fsetXbInv;
      atlas::FieldSet fsetFgInv;
      for (const auto & var : saberOuterBlockParams.inverseVars.value().variables()) {
        fsetXbInv.add(fsetXb.field(var));
        fsetFgInv.add(fsetFg.field(var));
      }

      // Apply left inverse
      saberBlockChain.lastOuterBlock().leftInverseMultiply(fsetXbInv);
      saberBlockChain.lastOuterBlock().leftInverseMultiply(fsetFgInv);

      // Copy the fields back
      for (const auto & var : saberOuterBlockParams.inverseVars.value().variables()) {
        fsetXb.field(var) = fsetXbInv.field(var);
        fsetFg.field(var) = fsetFgInv.field(var);
      }
    }

    // Adjoint test
    if (covarConf.getBool("adjoint test")) {
      // Get tolerance
      const double localAdjointTolerance =
        saberOuterBlockParams.adjointTolerance.value().get_value_or(
        covarConf.getDouble("adjoint tolerance"));

      // Run test
      saberBlockChain.lastOuterBlock().adjointTest(geom.getComm(),
                                                   outerGeometryData.back().get(),
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
                          << saberBlockChain.lastOuterBlock().blockName() << std::endl;
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
        saberBlockChain.lastOuterBlock().inverseTest(innerGeometryData,
                                                     activeInnerVars,
                                                     outerGeometryData.back().get(),
                                                     activeOuterVars,
                                                     innerVarsToCompare,
                                                     outerVarsToCompare,
                                                     innerInverseTolerance,
                                                     outerInverseTolerance,
                                                     geom.timeComm().rank());
      }
    }

    // Update outer geometry and variables for the next block
    outerGeometryData.push_back(innerGeometryData);
    outerVars = innerVars;

    // Apply adjoint to progressively get central fields
    saberBlockChain.lastOuterBlock().multiplyAD(saberBlockChain.centralFieldSet());
  }

  oops::Log::trace() << "buildOuterBlocks done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void buildCentralBlock(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResolutionGeom,
                       const oops::GeometryData & outerGeometryData,
                       oops::Variables & outerVars,
                       const atlas::FieldSet & fsetXb,
                       const atlas::FieldSet & fsetFg,
                       const util::DateTime & validTimeOfXbFg,
                       std::vector<atlas::FieldSet> & fsetEns,
                       std::vector<atlas::FieldSet> & dualResolutionFsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const SaberCentralBlockParametersWrapper & saberCentralBlockParamsWrapper,
                       SaberBlockChain & saberBlockChain) {
  oops::Log::trace() << "buildCentralBlock starting" << std::endl;

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Get active variables
  oops::Variables activeVars = getActiveVars(saberCentralBlockParams, outerVars);

  // Check that active variables are present in variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(outerVars.has(var));
  }

  // Create central block
  saberBlockChain.centralBlockInit(SaberCentralBlockFactory::create(
                                   outerGeometryData,
                                   activeVars,
                                   covarConf,
                                   saberCentralBlockParams,
                                   fsetXb,
                                   fsetFg,
                                   validTimeOfXbFg,
                                   geom.timeComm().rank()));

  // Read and add model fields
  saberBlockChain.centralBlock().read(geom, outerVars, validTimeOfXbFg);

  if (saberCentralBlockParams.doCalibration()) {
    // Block calibration

    // Ensemble configuration
    eckit::LocalConfiguration ensembleConf
      = covarConf.getSubConfiguration("ensemble configuration");

    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      saberBlockChain.centralBlock().iterativeCalibrationInit();

      // Get ensemble size
      size_t nens = ensembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < nens; ++ie) {
        // Read ensemble member
        atlas::FieldSet fset;
        readEnsembleMember(geom,
                           saberBlockChain.incVars(),
                           validTimeOfXbFg,
                           ensembleConf,
                           ie,
                           fset);

        // Apply outer blocks inverse (all of them)
        oops::Log::info() << "Info     : Apply outer blocks inverse (all of them)" << std::endl;
        saberBlockChain.leftInverseMultiply(fset);

        // Use FieldSet in the central block
        oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
        saberBlockChain.centralBlock().iterativeCalibrationUpdate(fset);
      }

      // Finalization
      oops::Log::info() << "Info     : Finalization" << std::endl;
      saberBlockChain.centralBlock().iterativeCalibrationFinal();
    } else {
      // Ensemble transform and localization for the Ensemble central block
      if (saberCentralBlockParams.saberBlockName.value() == "Ensemble") {
        // Read inflation field
        oops::Log::info() << "Info     : Read inflation field" << std::endl;
        saberBlockChain.centralBlock().read();

        // Apply inflation on ensemble members
        oops::Log::info() << "Info     : Apply inflation on ensemble members" << std::endl;
        saberBlockChain.centralBlock().applyInflation(fsetEns);

        // Ensemble transform
        const auto & ensTransConf = saberCentralBlockParams.ensembleTransform.value();
        if (ensTransConf != boost::none) {
          // Initialize ensemble transform blockchain
          std::unique_ptr<SaberBlockChain> ensTransBlockChain;
          ensTransBlockChain = ensembleBlockChain<MODEL>(geom,
                                                         outerVars,
                                                         fsetXb,
                                                         fsetXb,
                                                         validTimeOfXbFg,
                                                         *ensTransConf);

          // Left inverse of ensemble transform on ensemble members
          oops::Log::info() << "Info     : Left inverse of ensemble transform on ensemble members"
                            << std::endl;
          for (auto & fset : fsetEns) {
            for (auto itob = ensTransBlockChain->outerBlocks().begin();
              itob != ensTransBlockChain->outerBlocks().end(); ++itob) {
                itob->get()->leftInverseMultiply(fset);
            }
          }

          // Add ensemble transform blocks to outer blocks
          oops::Log::info() << "Info     : Add ensemble transform blocks to outer blocks"
                            << std::endl;
          std::move(ensTransBlockChain->outerBlocks().begin(),
            ensTransBlockChain->outerBlocks().end(),
            std::back_inserter(saberBlockChain.outerBlocks()));
          ensTransBlockChain->outerBlocks().erase(ensTransBlockChain->outerBlocks().begin(),
            ensTransBlockChain->outerBlocks().end());
        }

        // Localization
        const auto & locConf = saberCentralBlockParams.localization.value();
        if (locConf != boost::none) {
          // Initialize localization blockchain
          std::unique_ptr<SaberBlockChain> locBlockChain;
          locBlockChain = ensembleBlockChain<MODEL>(geom,
                                                    outerVars,
                                                    fsetXb,
                                                    fsetXb,
                                                    validTimeOfXbFg,
                                                    *locConf);
          // Move localization blockchain
          oops::Log::info() << "Info     : Move localization blockchain" << std::endl;
          saberBlockChain.centralBlock().setLocalization(std::move(locBlockChain));
        }
      }

      // Direct calibration
      oops::Log::info() << "Info     : Direct calibration" << std::endl;
      saberBlockChain.centralBlock().directCalibration(fsetEns);
    }
  } else if (saberCentralBlockParams.doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    saberBlockChain.centralBlock().read();
  }

  // Dual resolution ensemble
  if (covarConf.has("dual resolution ensemble configuration")) {
    oops::Log::info() << "Info     : Dual resolution setup" << std::endl;

    // Dual resolution setup
    saberBlockChain.centralBlock().dualResolutionSetup(dualResolutionGeom.generic());

    // Ensemble configuration
    eckit::LocalConfiguration dualResolutionEnsembleConf
      = covarConf.getSubConfiguration("dual resolution ensemble configuration");

    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      saberBlockChain.centralBlock().iterativeCalibrationInit();

      // Get dual resolution ensemble size
      size_t dualResolutionNens = dualResolutionEnsembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < dualResolutionNens; ++ie) {
        // Read ensemble member
        atlas::FieldSet fset;
        readEnsembleMember(dualResolutionGeom,
                           saberBlockChain.incVars(),
                           validTimeOfXbFg,
                           dualResolutionEnsembleConf,
                           ie,
                           fset);

        // Use FieldSet in the central block
        oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
        saberBlockChain.centralBlock().iterativeCalibrationUpdate(fset);
     }

      // Finalization
      oops::Log::info() << "Info     : Finalization" << std::endl;
      saberBlockChain.centralBlock().iterativeCalibrationFinal();
    } else {
      // Direct calibration
      oops::Log::info() << "Info     : Direct calibration" << std::endl;
      saberBlockChain.centralBlock().directCalibration(dualResolutionFsetEns);
    }
  }

  // Write calibration data
  if (saberCentralBlockParams.doCalibration()) {
    oops::Log::info() << "Info     : Write calibration data" << std::endl;
    saberBlockChain.centralBlock().write(geom, outerVars, validTimeOfXbFg);
    saberBlockChain.centralBlock().write();
  }

  // Write final ensemble
  if (covarConf.has("output ensemble")) {
    // Write ensemble
    OutputEnsembleParameters<MODEL> outputEnsemble;
    outputEnsemble.deserialize(covarConf.getSubConfiguration("output ensemble"));
    writeEnsemble(geom,
                  activeVars,
                  validTimeOfXbFg,
                  outputEnsemble,
                  covarConf.getSubConfiguration("ensemble configuration"),
                  iterativeEnsembleLoading,
                  saberBlockChain,
                  fsetEns);
  }

  // Adjoint test
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    saberBlockChain.centralBlock().adjointTest(geom.getComm(),
                                               outerGeometryData,
                                               activeVars,
                                               localAdjointTolerance);
  }

  oops::Log::trace() << "buildCentralBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
