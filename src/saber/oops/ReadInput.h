/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void readInputFields(
  const oops::Geometry<MODEL> & geom,
  const oops::Variables & vars,
  const util::DateTime & date,
  const std::vector<eckit::LocalConfiguration> & inputFieldConfs,
  std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << "readInputFields starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(geom.getComm().size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Read block input fields function
  if (inputFieldConfs.size() > 0) {
    // Loop over block input fields
    for (const auto & inputFieldConf : inputFieldConfs) {
      // Get input field file configuration
      eckit::LocalConfiguration file = inputFieldConf.getSubConfiguration("file");

      // Replace patterns
      util::seekAndReplace(file, "_MPI_", mpi);
      util::seekAndReplace(file, "_OMP_", omp);

      // Read block input field as Increment
      oops::Increment<MODEL> dx(geom, vars, date);
      dx.read(file);

      // Define FieldSet name
      std::string name = inputFieldConf.getString("parameter");
      if (inputFieldConf.has("component")) {
        name += "::" + std::to_string(inputFieldConf.getInt("component"));
      }

      // Transform Increment into FieldSet
      oops::Log::test() << "Norm of input parameter " << name << ": " << dx.norm() << std::endl;
      atlas::FieldSet fset;
      fset.name() = name;
      for (const auto field : dx.fieldSet()) {
        fset.add(field);
      }
      fsetVec.push_back(fset);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void readEnsemble(
  const oops::Geometry<MODEL> & geom,
  const oops::Variables & vars,
  const oops::State<MODEL> & xb,
  const oops::State<MODEL> & fg,
  const SaberBlockParametersBase & params,
  std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << "readEnsemble starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(geom.getComm().size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Ensemble pointer
  std::unique_ptr<oops::IncrementEnsemble<MODEL>> ensemble_;

  // Ensemble of states, perturbation using the mean
  eckit::LocalConfiguration ensembleConf = params.ensemble.value()
    .get_value_or(eckit::LocalConfiguration());
  if (!ensembleConf.empty()) {
    oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                      << std::endl;
    util::seekAndReplace(ensembleConf, "_MPI_", mpi);
    util::seekAndReplace(ensembleConf, "_OMP_", omp);
    oops::IncrementEnsembleFromStatesParameters<MODEL> ensembleParams;
    ensembleParams.validateAndDeserialize(ensembleConf);
    ensemble_.reset(new oops::IncrementEnsemble<MODEL>(ensembleParams, xb, fg, geom, vars));
  }

  // Increment ensemble from increments on disk
  eckit::LocalConfiguration ensemblePertConf = params.ensemblePert.value()
    .get_value_or(eckit::LocalConfiguration());
  if (!ensemblePertConf.empty()) {
    util::seekAndReplace(ensemblePertConf, "_MPI_", mpi);
    util::seekAndReplace(ensemblePertConf, "_OMP_", omp);
    oops::IncrementEnsembleParameters<MODEL> ensemblePertParams;
    ensemblePertParams.validateAndDeserialize(ensemblePertConf);
    oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
    ensemble_.reset(new oops::IncrementEnsemble<MODEL>(geom, vars, ensemblePertParams));
  }

  // Increment ensemble from difference of two states
  eckit::LocalConfiguration ensembleBaseConf = params.ensembleBase.value()
    .get_value_or(eckit::LocalConfiguration());
  eckit::LocalConfiguration ensemblePairsConf = params.ensemblePairs.value()
    .get_value_or(eckit::LocalConfiguration());
  if (!ensembleBaseConf.empty() && !ensemblePairsConf.empty()) {
    oops::Log::info() << "Info     : Increment ensemble from difference of two states"
                      << std::endl;
    util::seekAndReplace(ensembleBaseConf, "_MPI_", mpi);
    util::seekAndReplace(ensembleBaseConf, "_OMP_", omp);
    oops::StateEnsembleParameters<MODEL> ensembleBaseParams;
    ensembleBaseParams.validateAndDeserialize(ensembleBaseConf);
    util::seekAndReplace(ensemblePairsConf, "_MPI_", mpi);
    util::seekAndReplace(ensemblePairsConf, "_OMP_", omp);
    oops::StateEnsembleParameters<MODEL> ensemblePairsParams;
    ensemblePairsParams.validateAndDeserialize(ensemblePairsConf);
    ensemble_.reset(new oops::IncrementEnsemble<MODEL>(geom, vars, ensembleBaseParams,
      ensemblePairsParams));
  }

  if (ensemble_) {
    for (unsigned int ie = 0; ie < ensemble_->size(); ++ie) {
       // Potentially apply some transforms on ensemble members
       // TODO(Benjamin)

       // Transform Increment into FieldSet
       atlas::FieldSet fset;
       fset.name() = "ensemble member";
       for (const auto field : (*ensemble_)[ie].fieldSet()) {
         fset.add(field);
       }
      fsetVec.push_back(fset);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
