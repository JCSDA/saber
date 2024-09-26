/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP.h"

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/FieldSets.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/bump/type_bump_parameters.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

BUMP::BUMP(const oops::GeometryData & geometryData,
           const oops::Variables & vars,
           const eckit::Configuration & covarConf,
           const BUMPParameters & params,
           const eckit::LocalConfiguration & fieldsMetaData,
           const oops::FieldSet3D & xb)
  : keyBUMP_(), comm_(geometryData.comm()), fspace_(geometryData.functionSpace()), vars_(vars),
  validTime_(xb.validTime()), covarConf_(covarConf), bumpConf_(params.toConfiguration()),
  nens_(), waitForDualResolution_(false), gridUid_(util::getGridUid(geometryData.functionSpace())),
  dualResolutionGridUid_("") {
  oops::Log::trace() << classname() << "::BUMP starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm_.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  oops::Log::info() << "Info     : +++ MPI tasks:      " << mpi << std::endl;
  oops::Log::info() << "Info     : +++ OpenMP threads: " << omp << std::endl;

  // Initialization
  nens_.reserve(2);
  nens_[0] = 0;
  if (covarConf_.has("ensemble configuration")) {
    nens_[0] = covarConf_.getSubConfiguration("ensemble configuration").getInt("ensemble size");
  }
  nens_[1] = 0;
  if (covarConf_.has("dual resolution ensemble configuration")) {
    nens_[1] = covarConf_.getSubConfiguration("dual resolution ensemble configuration")
      .getInt("ensemble size");
    waitForDualResolution_ = true;
  }
  iterativeEnsembleLoading_ = covarConf_.getBool("iterative ensemble loading", false);

  // Case where size are specified in the BUMP configuration
  // TODO(Benjamin): when is this necessary?
  if ((nens_[0] == 0) && bumpConf_.has("ensemble sizes.total ensemble size")) {
    // Ensemble 1 size from configuration
    nens_[0] = bumpConf_.getInt("ensemble sizes.total ensemble size");
  }
  if ((nens_[1] == 0) && bumpConf_.has("ensemble sizes.total lowres ensemble size")) {
    // Ensemble 1 size from configuration
    nens_[1] = bumpConf_.getInt("ensemble sizes.total lowres ensemble size");
  }

  // Update ensemble sizes
  bumpConf_.set("ensemble sizes.total ensemble size", nens_[0]);
  bumpConf_.set("ensemble sizes.total lowres ensemble size", nens_[1]);

  // Set iterative ensemble loading flag
  bumpConf_.set("external.iterative algorithm", iterativeEnsembleLoading_);

  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get the grids configuration from input configuration and complete it
  if (bumpConf_.has("grids")) {
    // Get grids from input configuration
    bumpConf_.get("grids", grids);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Check grids number
  if (grids.size() == 0) {
    oops::Log::info() << "BUMP: grid size is zero" << std::endl;
    std::abort();
  }

  // Loop over grids
  size_t jgrid = 0;
  for (auto & grid : grids) {
    // Merge bumpConf_ into grid
    grid = util::mergeConfigs(grid, bumpConf_);

    // Replace patterns
    util::seekAndReplace(grid, "_MPI_", mpi);
    util::seekAndReplace(grid, "_OMP_", omp);

    // New files
    if (jgrid > 0) {
      grid.set("io.new files", false);
    }

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grid.has("model.variables")) {
      grid.get("model.variables", vars_str);
    }
    if (vars_str.size() == 0) {
      vars_str = vars_.variables();
      grid.set("model.variables", vars_str);
    }

    // Get the number of levels and the 2D variables
    int nl0 = 0;
    std::vector<std::string> var2d;
    for (const auto & var : vars_str) {
      bool varFound = false;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if (var == vars_[jvar].name()) {
          varFound = true;
          int nl0_tmp = static_cast<int>(vars_[var].getLevels());
          if (nl0 > 1) {
            // Check that nl0_tmp is either 1 or nl0
            if ((nl0_tmp != 1) && (nl0_tmp != nl0)) {
             oops::Log::info() << "BUMP::BUMP: inconsistent number of levels in BUMP" << std::endl;
              std::abort();
            }
          }
          nl0 = std::max(nl0, nl0_tmp);

          // 2D variable flag
          if (nl0_tmp == 1) {
            var2d.push_back(var);
          }
        }
      }
      if (!varFound) {
        oops::Log::info() << "BUMP: inconsistent variable names" << std::endl;
        std::abort();
      }
    }
    grid.set("model.nl0", nl0);
    grid.set("model.2d variables", var2d);

    // Add level index for 2D fields
    if (!grid.has("model.level for 2d variables")) {
      ModelDef def;
      grid.set("model.level for 2d variables", def.lev2d.second);
    }

    // Add vertical coordinate name
    std::string vertCoordName;
    for (const auto & var : vars_str) {
      if (var2d.size() == vars_str.size() || vars_[var].getLevels() > 1) {
        const std::string key = var + ".vert_coord";
        if (vertCoordName.empty()) {
          vertCoordName = fieldsMetaData.getString(key, "");
        } else {
          ASSERT(fieldsMetaData.getString(key, "vert_coord") == vertCoordName);
        }
      }
    }
    if (vertCoordName.empty() && geometryData.fieldSet().has("vert_coord")) {
      vertCoordName = "vert_coord";
    }
    grid.set("external.vertical coordinate name", vertCoordName);

    // Add geographical mask name
    std::string gmaskName;
    for (const auto & var : vars_str) {
      if (var2d.size() == vars_str.size() || vars_[var].getLevels() > 1) {
        const std::string key = var + ".gmask";
        if (gmaskName.empty()) {
          gmaskName = fieldsMetaData.getString(key, "");
        } else {
          ASSERT(fieldsMetaData.getString(key, "gmask") == gmaskName);
        }
      }
    }
    if (gmaskName.empty() && geometryData.fieldSet().has("gmask")) {
      gmaskName = "gmask";
    }
    grid.set("external.geographical mask name", gmaskName);

    // Create BUMP instance
    oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    oops::Log::info() << "Info     : +++ Create BUMP instance " << (jgrid+1) << " / "
      << grids.size()
                  << std::endl;
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, &comm_, fspace_.get(), geometryData.fieldSet().get(), grid,
      &oops::LibOOPS::instance().infoChannel(), &oops::LibOOPS::instance().testChannel());
    keyBUMP_.push_back(keyBUMP);
    ++jgrid;
  }

  // Get max number of components
  size_t ncmp = 0;
  if (bumpConf_.has("input atlas files")) {
    for (const auto & input : bumpConf_.getSubConfigurations("input atlas files")) {
      size_t icmp = input.getInt("component", 1);
      ncmp = std::max(ncmp, icmp);
    }
  }
  if (bumpConf_.has("input model files")) {
    for (const auto & input : bumpConf_.getSubConfigurations("input model files")) {
      size_t icmp = input.getInt("component", 1);
      ncmp = std::max(ncmp, icmp);
    }
  }

  if (ncmp > 0) {
    // Set number of components
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      oops::Log::info() << "Info     : +++ Set number of components for BUMP instance " << (jgrid+1)
                    << " / " << keyBUMP_.size() << std::endl;
      bump_set_ncmp_f90(keyBUMP_[jgrid], ncmp);
    }
  }
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  oops::Log::info() << "Info     : +++ End of BUMP creation" << std::endl;
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;

  oops::Log::trace() << classname() << "::BUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP::~BUMP() {
  oops::Log::trace() << classname() << "::~BUMP starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (keyBUMP_[jgrid] > 0) {
      oops::Log::info() << "Info     :"
                    << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      oops::Log::info() << "Info     : +++ Terminate BUMP instance " << (jgrid+1) << " / "
                    << keyBUMP_.size() << std::endl;
      bump_dealloc_f90(keyBUMP_[jgrid]);
    }
  }
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  oops::Log::info() << "Info     : +++ All BUMP instances terminated" << std::endl;
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;

  oops::Log::trace() << classname() << "::~BUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::readAtlasFiles() {
  oops::Log::trace() << classname() << "::readAtlasFiles starting" << std::endl;

  if (bumpConf_.has("input atlas files")) {
    std::vector<eckit::LocalConfiguration> inputAtlasFilesConf
      = bumpConf_.getSubConfigurations("input atlas files");
    for (const auto & conf : inputAtlasFilesConf) {
      // Get file configuration
      eckit::LocalConfiguration file = this->getFileConf(comm_, conf);

      // Get parameter and component
      const std::string param = conf.getString("parameter");
      const int icmp = conf.getInt("component", 1);

      // Create FieldSet
      oops::FieldSet3D fset(validTime_, comm_);
      fset.name() = param + " - " + std::to_string(icmp);

      // Read file
      fset.read(fspace_, vars_, file);

      // Print FieldSet norm
      oops::Log::test() << "+++ Norm of input parameter " << fset.name() << ": "
                        << fset.norm(vars_)
                        << std::endl;

      // Add fields
      this->addField(fset);
    }
  }

  oops::Log::trace() << classname() << "::readAtlasFiles done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> BUMP::getReadConfs(
  const std::vector<eckit::LocalConfiguration> confs) {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  for (const auto & conf : confs) {
    // Get parameter and component
    const std::string param = conf.getString("parameter");
    const int icmp = conf.getInt("component", 1);
    const std::string name = param + " - " + std::to_string(icmp);

    // Get file configuration
    eckit::LocalConfiguration file = this->getFileConf(comm_, conf);

    // Add pair
    inputs.push_back(std::make_pair(name, file));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void BUMP::addField(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::addField starting" << std::endl;

  // Check fset grid UID
  if (fset.size() > 0) {
    if (fset.getGridUid() != gridUid_) {
      oops::Log::info() << "BUMP: wrong grid UID" << std::endl;
      std::abort();
    }
  }

  // Set parameter
  size_t pos = fset.name().find(" - ");
  size_t icmp = 1;
  std::string param = fset.name();
  if (pos != std::string::npos) {
    icmp = std::stoi(fset.name().substr(pos+3));
    param = fset.name().substr(0, pos);
  }
  const int npar = param.size();
  const char *cpar = param.c_str();

  // Set parameter
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_set_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, fset.get());
  }

  oops::Log::trace() << classname() << "::addField done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::addEnsemble(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::addEnsemble starting" << std::endl;

  // Initialize ensemble index
  size_t ie = 0;
  for (size_t jj = 0; jj < fsetEns.ens_size(); ++jj) {
    // Get geometry index (iterative update should be done after for the dual resolution part)
    // and check grid UID
    size_t igeom;
    if (dualResolutionGridUid_ == "") {
      igeom = 0;
      if (fsetEns[jj].getGridUid() != gridUid_) {
        oops::Log::info() << "BUMP::addEnsemble: wrong grid UID" << std::endl;
        std::abort();
      }
    } else {
      igeom = 1;
      if (fsetEns[jj].getGridUid() != dualResolutionGridUid_) {
        oops::Log::info() << "BUMP::addEnsemble: wrong dual resolution grid UID" << std::endl;
        std::abort();
      }
    }

    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      // Initial message
      if (ie == 0) {
        oops::Log::info() << "Info     :"
                      << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      << std::endl;
        oops::Log::info() << "Info     : +++ Add members of ensemble " << (igeom+1) << std::endl;
      }

      // Add member
      oops::Log::info() << "Info     :        Member " << ie+1 << " / " << nens_[igeom]
        << std::endl;
      bump_add_member_f90(keyBUMP_[jgrid], fsetEns[jj].fieldSet().get(), ie+1, igeom+1);
    }

    // Update member index
    ++ie;
  }

  oops::Log::trace() << classname() << "::addEnsemble starting" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::dualResolutionSetup(const atlas::FunctionSpace & fspace,
                               const atlas::FieldSet & fields) {
  oops::Log::trace() << classname() << "::dualResolutionSetup starting" << std::endl;

  // Set dual resolution grid UID
  dualResolutionGridUid_ = util::getGridUid(fspace);

  // Dual resolution setup
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_dual_resolution_setup_f90(keyBUMP_[jgrid], fspace.get(), fields.get());
  }

  oops::Log::trace() << classname() << "::dualResolutionSetup done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::iterativeUpdate(const oops::FieldSet3D & fset,
                           const size_t & ie) {
  oops::Log::trace() << classname() << "::iterativeUpdate starting" << std::endl;

  // Get geometry index (iterative update should be done after for the dual resolution part)
  // and check grid UID
  size_t igeom;
  if (dualResolutionGridUid_ == "") {
    igeom = 0;
    if (fset.getGridUid() != gridUid_) {
      oops::Log::info() << "BUMP::iterativeUpdate: wrong grid UID" << std::endl;
      std::abort();
    }
  } else {
    igeom = 1;
    if (fset.getGridUid() != dualResolutionGridUid_) {
      oops::Log::info() << "BUMP::iterativeUpdate: wrong dual resolution grid UID" << std::endl;
      std::abort();
    }
  }

  // Print info
  oops::Log::info() << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  oops::Log::info() << "Info     : +++ Load member " << ie+1 << " / " << nens_[igeom] << std::endl;

  // Get driver keys
  bool new_vbal_cov = bumpConf_.getBool("drivers.compute vertical covariance", false);
  bool new_var = bumpConf_.getBool("drivers.compute variance", false);
  bool new_mom = bumpConf_.getBool("drivers.compute moments", false);
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (igeom == 0) {
      if (new_vbal_cov) {
        // Update vertical covariance
        bump_update_vbal_cov_f90(keyBUMP_[jgrid], fset.get(), ie+1);
      }
      if (new_var) {
        // Update variance
        bump_update_var_f90(keyBUMP_[jgrid], fset.get(), ie+1);
      }
    }
    if (new_mom) {
      // Update moments
      bump_update_mom_f90(keyBUMP_[jgrid], fset.get(), ie+1, igeom+1);
    }
  }

  // Print info
  if (ie+1 == nens_[igeom]) {
    oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    oops::Log::info() << "Info     : +++ End of iterative update" << std::endl;
    oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
  }

  oops::Log::trace() << classname() << "::iterativeUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::writeAtlasFiles() const {
  oops::Log::trace() << classname() << "::writeAtlasFiles starting" << std::endl;

  if (bumpConf_.has("output atlas files")) {
    std::vector<eckit::LocalConfiguration> outputAtlasFilesConf
      = bumpConf_.getSubConfigurations("output atlas files");
    for (const auto & output : this->fieldsToWrite(outputAtlasFilesConf)) {
      // Get file configuration
      eckit::LocalConfiguration file = this->getFileConf(comm_, output.first);

      // Print FieldSet norm
      oops::Log::test() << "+++ Norm of output parameter " << output.second.name() << ": "
                        << output.second.norm(vars_)
                        << std::endl;

      // Write FieldSet
      output.second.write(file);
    }
  }

  oops::Log::trace() << classname() << "::writeAtlasFiles done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> BUMP::fieldsToWrite(
  const std::vector<eckit::LocalConfiguration> confs) const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Define outputs vector
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs;

  // Write output fields
  for (const auto & conf : confs) {
    // Get parameter and component
    const std::string param = conf.getString("parameter");
    const int icmp = conf.getInt("component", 1);

    // Create FieldSet
    oops::FieldSet3D fset(validTime_, comm_);
    fset.name() = param + " - " + std::to_string(icmp);

    // Get FieldSet
    const int npar = param.size();
    const char *cpar = param.c_str();
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      bump_get_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, fset.get());
    }

    // Get file configuration
    eckit::LocalConfiguration file = this->getFileConf(comm_, conf);

    // Add pair
    outputs.push_back(std::make_pair(file, fset));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void BUMP::runDrivers() {
  oops::Log::trace() << classname() << "::runDrivers starting" << std::endl;

  if (waitForDualResolution_) {
    waitForDualResolution_ = false;
  } else {
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      oops::Log::info() << "Info     :"
                    << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      oops::Log::info() << "Info     : +++ Run drivers for BUMP instance " << (jgrid+1)
                    << " / " << keyBUMP_.size() << std::endl;
      bump_run_drivers_f90(keyBUMP_[jgrid]);
      bump_partial_dealloc_f90(keyBUMP_[jgrid]);
    }
    oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    oops::Log::info() << "Info     : +++ End of BUMP drivers" << std::endl;
    oops::Log::info() << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
  }

  oops::Log::trace() << classname() << "::runDrivers done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbal(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyVbal starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::done starting" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbalAd(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyVbalAd starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::multiplyVbalAd done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyVbal(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyVbal starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::inverseMultiplyVbal done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyStdDev(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyStdDev starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::multiplyStdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyStdDev(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyStdDev starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::inverseMultiplyStdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::randomizeNicas(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomizeNicas starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_randomize_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::randomizeNicas done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicas(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyNicas starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::multiplyNicas done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUV(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyPsiChiToUV starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::multiplyPsiChiToUV done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUVAd(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyPsiChiToUVAd starting" << std::endl;

  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_ad_f90(keyBUMP_[jgrid], fset.get());
  }

  oops::Log::trace() << classname() << "::multiplyPsiChiToUVAd done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t BUMP::getCvSize() const {
  oops::Log::trace() << classname() << "::getCvSize starting" << std::endl;

  size_t n = 0;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    n += nTmp;
  }

  oops::Log::trace() << classname() << "::getCvSize done" << std::endl;
  return n;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicasSqrt(const atlas::Field & cv,
                             oops::FieldSet3D & fset,
                             const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplyNicasSqrt starting" << std::endl;

  int index = offset;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_sqrt_f90(keyBUMP_[jgrid], cv.get(), fset.get(), index);
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    index += nTmp;
  }

  oops::Log::trace() << classname() << "::multiplyNicasSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicasSqrtAd(const oops::FieldSet3D & fset,
                               atlas::Field & cv,
                               const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplyNicasSqrtAd starting" << std::endl;

  int index = offset;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_sqrt_ad_f90(keyBUMP_[jgrid], fset.get(), cv.get(), index);
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    index += nTmp;
  }

  oops::Log::trace() << classname() << "::multiplyNicasSqrtAd done" << std::endl;
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration BUMP::getFileConf(const eckit::mpi::Comm & comm,
                                            const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::getFileConf starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm_.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Get IO configuration
  eckit::LocalConfiguration file = conf.getSubConfiguration("file");

  // Replace patterns
  util::seekAndReplace(file, "_MPI_", mpi);
  util::seekAndReplace(file, "_OMP_", omp);

  oops::Log::trace() << classname() << "::getFileConf done" << std::endl;
  return file;
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
