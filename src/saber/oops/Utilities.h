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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/ErrorCovarianceParameters.h"

namespace saber {

oops::Variables getActiveVars(const SaberBlockParametersBase &,
                              const oops::Variables &);

void setMember(eckit::LocalConfiguration & conf,
               const int & member);

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi);

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

}  // namespace saber
