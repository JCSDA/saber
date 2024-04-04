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

#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/StateSet.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelData.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/ErrorCovarianceParameters.h"

namespace oops {
  class FieldSet3D;
}

namespace saber {

// -----------------------------------------------------------------------------

oops::Variables getActiveVars(const SaberBlockParametersBase & params,
                              const oops::Variables & defaultVars);

// -----------------------------------------------------------------------------

oops::Variables getUnionOfInnerActiveAndOuterVars(const SaberBlockParametersBase & params,
                                                  const oops::Variables & outerVars);

// -----------------------------------------------------------------------------

oops::Variables getInnerOnlyVars(const SaberBlockParametersBase & params,
                                 const oops::Variables & outerVars);

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

// -----------------------------------------------------------------------------

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::Variables & vars);

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::Variables & varsToAllocate,
                           const oops::Variables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace);

// -----------------------------------------------------------------------------

template<typename MODEL>
oops::FieldSets readEnsemble(const oops::Geometry<MODEL> & geom,
                                      const oops::Variables & modelvars,
                                      const oops::State4D<MODEL> & xb,
                                      const oops::State4D<MODEL> & fg,
                                      const eckit::LocalConfiguration & inputConf,
                                      const bool & iterativeEnsembleLoading,
                                      eckit::LocalConfiguration & outputConf) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Prepare ensemble configuration
  oops::Log::info() << "Info     : Prepare ensemble configuration" << std::endl;

  // Fill output configuration and set ensemble size
  size_t nens = 0;
  size_t ensembleFound = 0;
  eckit::LocalConfiguration varConf;

  // Ensemble of states, perturbation using the mean
  oops::IncrementEnsembleFromStatesParameters<MODEL> ensembleParams;
  eckit::LocalConfiguration ensembleConf;
  if (inputConf.has("ensemble")) {
    ensembleConf = inputConf.getSubConfiguration("ensemble");
    ensembleParams.deserialize(ensembleConf);
    nens = ensembleParams.states.size();
    outputConf.set("ensemble", ensembleConf);
    varConf = ensembleParams.states.getStateConfig(0, 0);
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
    ensemblePertParams.getIncrementParameters(0).serialize(varConf);
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
    varConf = ensembleBaseParams.getStateConfig(0, 0);
    outputConf.set("ensemble base", ensembleBase);
    outputConf.set("ensemble pairs", ensemblePairs);
    ++ensembleFound;
  }

  // Set ensemble size
  outputConf.set("ensemble size", nens);

  // Check number of ensembles in yaml
  ASSERT(ensembleFound <= 1);

  oops::Variables vars(varConf.has("variables") ?
    oops::Variables{varConf.getStringVector("variables")} :
    modelvars);

  if (!iterativeEnsembleLoading) {
    // Full ensemble loading
    oops::Log::info() << "Info     : Read full ensemble" << std::endl;

    // Ensemble of states, perturbation using the mean
    if (!ensembleConf.empty()) {
      oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                        << std::endl;
      oops::StateSet<MODEL> tmp(geom, ensembleParams.toConfiguration(), xb.commTime());
      oops::IncrementSet<MODEL> ensemble(geom, vars, tmp, true);
      ensemble -= ensemble.ens_mean();
      oops::FieldSets fsetEns(ensemble);
      return fsetEns;
    }

    // Increment ensemble from increments on disk
    if (!ensemblePert.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
      oops::IncrementSet<MODEL> ensemble(geom, vars, xb.times(),
                                         ensemblePertParams.toConfiguration(), xb.commTime());
      oops::FieldSets fsetEns(ensemble);
      return fsetEns;
    }

    // Increment ensemble from difference of two states
    if (!ensembleBase.empty() && !ensemblePairs.empty()) {
      oops::Log::info() << "Info     : Increment ensemble from difference of two states"
                        << std::endl;
      oops::StateSet<MODEL> states1(geom, ensembleBaseParams.toConfiguration(), xb.commTime());
      oops::StateSet<MODEL> states2(geom, ensemblePairsParams.toConfiguration(), xb.commTime());
      oops::IncrementSet<MODEL> ensemble(geom, vars, states1.times(), states1.commTime(),
                                         states1.members(), states1.commEns());
      ensemble.diff(states1, states2);
      oops::FieldSets fsetEns(ensemble);
      return fsetEns;
    }
  }
  // Return empty ensemble if none was returned before
  std::vector<util::DateTime> dates;
  std::vector<int> ensmems;
  oops::FieldSets ensemble(dates, xb.commTime(), ensmems, xb.commEns());
  return ensemble;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readHybridWeight(const oops::Geometry<MODEL> & geom,
                      const oops::Variables & vars,
                      const util::DateTime & date,
                      const eckit::LocalConfiguration & conf,
                      oops::FieldSet3D & fset) {
  oops::Log::trace() << "readHybridWeight starting" << std::endl;

  oops::Log::info() << "Info     : Read hybrid weight" << std::endl;

  // Local copy
  eckit::LocalConfiguration localConf(conf);

  // Create Increment
  oops::Increment<MODEL> dx(geom, vars, date);

  // Read file
  dx.read(localConf);

  // Get FieldSet
  fset.shallowCopy(dx.fieldSet());

  oops::Log::trace() << "readHybridWeight done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readEnsembleMember(const oops::Geometry<MODEL> & geom,
                        const oops::Variables & vars,
                        const eckit::LocalConfiguration & conf,
                        const size_t & ie,
                        oops::FieldSet3D & fset) {
  oops::Log::trace() << "readEnsembleMember starting" << std::endl;

  oops::Log::info() << "Info     : Read ensemble member " << ie << std::endl;

  // Fill FieldSet
  size_t ensembleFound = 0;

  const size_t myrank = geom.timeComm().rank();

  if (conf.has("ensemble")) {
    // Ensemble of states passed as increments
    oops::StateEnsembleParameters<MODEL> states;
    states.deserialize(conf.getSubConfiguration("ensemble"));

    // Read state
    oops::State<MODEL> xx(geom, states.getStateConfig(ie, myrank));

    // Copy FieldSet
    fset.deepCopy(xx.fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble pert")) {
    // Increment ensemble from increments on disk
    oops::IncrementEnsembleParameters<MODEL> ensemblePertParams;
    ensemblePertParams.deserialize(conf.getSubConfiguration("ensemble pert"));

    // Read Increment
    oops::Increment<MODEL> dx(geom, vars, fset.validTime());
    dx.read(ensemblePertParams.getIncrementParameters(ie));

    // Get FieldSet
    fset.deepCopy(dx.fieldSet());

    ++ensembleFound;
  }

  if (conf.has("ensemble base") && conf.has("ensemble pairs")) {
    // Increment ensemble from difference of two states
    oops::StateEnsembleParameters<MODEL> ensembleBaseParams;
    ensembleBaseParams.deserialize(conf.getSubConfiguration("ensemble base"));
    oops::StateEnsembleParameters<MODEL> ensemblePairsParams;
    ensemblePairsParams.deserialize(conf.getSubConfiguration("ensemble pairs"));

    // Read states
    oops::State<MODEL> xxBase(geom, ensembleBaseParams.getStateConfig(ie, myrank));
    oops::State<MODEL> xxPairs(geom, ensemblePairsParams.getStateConfig(ie, myrank));

    // Compute difference
    oops::Increment<MODEL> dx(geom, vars, fset.validTime());
    dx.diff(xxPairs, xxBase);

    // Get FieldSet
    fset.deepCopy(dx.fieldSet());

    ++ensembleFound;
  }

  // Check number of ensembles in configuration
  ASSERT(ensembleFound <= 1);

  oops::Log::trace() << "readEnsembleMember done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
