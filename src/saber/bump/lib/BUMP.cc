/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/lib/BUMP.h"

#include <omp.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "saber/bump/lib/type_bump_parameters.h"


namespace bump_lib {

// -----------------------------------------------------------------------------

BUMP::BUMP(const eckit::mpi::Comm & comm,
           eckit::Channel & infoChannel,
           eckit::Channel & testChannel,
           const atlas::FunctionSpace & fspace,
           const atlas::FieldSet & fields,
           const std::vector<size_t> & variableSizes,
           const std::vector<std::string> & vars,
           const eckit::Configuration & covarConf,
           const eckit::Configuration & bumpConf) :
  keyBUMP_(), comm_(&comm), infoChannel_(&infoChannel), testChannel_(&testChannel),
  fspace_(fspace), variableSizes_(variableSizes), vars_(vars),
  covarConf_(covarConf), bumpConf_(bumpConf), nens_(), waitForDualResolution_(false),
  gridUid_(util::getGridUid(fspace)), dualResolutionGridUid_("") {
  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  *infoChannel_ << "Info     : +++ MPI tasks:      " << mpi << std::endl;
  *infoChannel_ << "Info     : +++ OpenMP threads: " << omp << std::endl;

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
    *infoChannel_ << "BUMP: grid size is zero" << std::endl;
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
      vars_str = vars_;
      grid.set("model.variables", vars_str);
    }

    // Get the number of levels and the 2D variables
    int nl0 = 0;
    std::vector<std::string> var2d;
    for (const auto & var : vars_str) {
      bool varFound = false;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if (var == vars_[jvar]) {
          varFound = true;
          int nl0_tmp = static_cast<int>(variableSizes_[jvar]);
          if (nl0 > 1) {
            // Check that nl0_tmp is either 1 or nl0
            if ((nl0_tmp != 1) && (nl0_tmp != nl0)) {
             *infoChannel_ << "BUMP::BUMP: inconsistent number of levels in BUMP" << std::endl;
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
        *infoChannel_ << "BUMP: inconsistent variable names" << std::endl;
        std::abort();
      }
    }
    grid.set("model.nl0", nl0);
    grid.set("model.2d variables", var2d);

    // Add level index for 2D fields
    if (!grid.has("model.lev2d")) {
      ModelDef def;
      grid.set("model.lev2d", def.lev2d.second);
    }

    // Create BUMP instance
    *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    *infoChannel_ << "Info     : +++ Create BUMP instance " << (jgrid+1) << " / " << grids.size()
                  << std::endl;
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, comm_, fspace_.get(), fields.get(),
                    grid, infoChannel_, testChannel_);
    keyBUMP_.push_back(keyBUMP);
    ++jgrid;
  }

  // Get max number of components
  size_t ncmp = 0;
  if (bumpConf.has("input atlas files")) {
    for (const auto & input : bumpConf.getSubConfigurations("input atlas files")) {
      size_t icmp = input.getInt("component", 1);
      ncmp = std::max(ncmp, icmp);
    }
  }
  if (bumpConf.has("input model files")) {
    for (const auto & input : bumpConf.getSubConfigurations("input model files")) {
      size_t icmp = input.getInt("component", 1);
      ncmp = std::max(ncmp, icmp);
    }
  }

  if (ncmp > 0) {
    // Set number of components
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      *infoChannel_ << "Info     : +++ Set number of components for BUMP instance " << (jgrid+1)
                    << " / " << keyBUMP_.size() << std::endl;
      bump_set_ncmp_f90(keyBUMP_[jgrid], ncmp);
    }
  }
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  *infoChannel_ << "Info     : +++ End of BUMP creation" << std::endl;
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
}

// -----------------------------------------------------------------------------

BUMP::~BUMP() {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (keyBUMP_[jgrid] > 0) {
      *infoChannel_ << "Info     :"
                    << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      *infoChannel_ << "Info     : +++ Terminate BUMP instance " << (jgrid+1) << " / "
                    << keyBUMP_.size() << std::endl;
      bump_dealloc_f90(keyBUMP_[jgrid]);
    }
  }
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  *infoChannel_ << "Info     : +++ All BUMP instances terminated" << std::endl;
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
}

// -----------------------------------------------------------------------------

void BUMP::readAtlasFiles() {
  if (bumpConf_.has("input atlas files")) {
    std::vector<eckit::LocalConfiguration> inputAtlasFilesConf
      = bumpConf_.getSubConfigurations("input atlas files");
    for (const auto & conf : inputAtlasFilesConf) {
      // Get file configuration
      eckit::LocalConfiguration file = this->getFileConf(*comm_, conf);

      // Get parameter and component
      const std::string param = conf.getString("parameter");
      const int icmp = conf.getInt("component", 1);

      // Create FieldSet
      atlas::FieldSet fset;
      fset.name() = param + " - " + std::to_string(icmp);

      // Read file
      util::readFieldSet(*comm_, fspace_, variableSizes_, vars_, file, fset);

      // Print FieldSet norm
        *testChannel_ << "+++ Norm of input parameter " << fset.name() << ": "
                      << util::normFieldSet(fset, vars_, *comm_)
                      << std::endl;

      // Add fields
      this->addField(fset);
    }
  }
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> BUMP::fieldsToRead(
  const std::vector<eckit::LocalConfiguration> confs) {
  // Write output fields
  for (const auto & conf : confs) {
    // Get file configuration
    eckit::LocalConfiguration file = this->getFileConf(*comm_, conf);

    // Get parameter and component
    const std::string param = conf.getString("parameter");
    const int icmp = conf.getInt("component", 1);

    // Create FieldSet
    atlas::FieldSet fset;
    fset.name() = param + " - " + std::to_string(icmp);

    // Add pair
    inputs_.push_back(std::make_pair(file, fset));
  }

  return inputs_;
}

// -----------------------------------------------------------------------------

void BUMP::addField(const atlas::FieldSet & fset) {
  // Check fset grid UID
  if (fset.size() > 0) {
    if (util::getGridUid(fset) != gridUid_) {
      *infoChannel_ << "BUMP: wrong grid UID" << std::endl;
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
}

// -----------------------------------------------------------------------------

void BUMP::addEnsemble(const std::vector<atlas::FieldSet> & fsetEns) {
  // Initialize ensemble index
  size_t ie = 0;

  for (const auto & fset : fsetEns) {
    // Get geometry index (iterative update should be done after for the dual resolution part)
    // and check grid UID
    size_t igeom;
    if (dualResolutionGridUid_ == "") {
      igeom = 0;
      if (util::getGridUid(fset) != gridUid_) {
        *infoChannel_ << "BUMP::iterativeUpdate: wrong grid UID" << std::endl;
        std::abort();
      }
    } else {
      igeom = 1;
      if (util::getGridUid(fset) != dualResolutionGridUid_) {
        *infoChannel_ << "BUMP::iterativeUpdate: wrong dual resolution grid UID" << std::endl;
        std::abort();
      }
    }

    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      // Initial message
      if (ie == 0) {
        *infoChannel_ << "Info     :"
                      << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      << std::endl;
        *infoChannel_ << "Info     : +++ Add members of ensemble " << (igeom+1) << std::endl;
      }

      // Add member
      *infoChannel_ << "Info     :        Member " << ie+1 << " / " << nens_[igeom] << std::endl;
      bump_add_member_f90(keyBUMP_[jgrid], fset.get(), ie+1, igeom+1);
    }

    // Update member index
    ++ie;
  }
}

// -----------------------------------------------------------------------------

void BUMP::dualResolutionSetup(const atlas::FunctionSpace & fspace,
                               const atlas::FieldSet & fields) {
  // Set dual resolution grid UID
  dualResolutionGridUid_ = util::getGridUid(fspace);

  // Dual resolution setup
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_dual_resolution_setup_f90(keyBUMP_[jgrid], fspace.get(), fields.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::iterativeUpdate(const atlas::FieldSet & fset, const size_t & ie) {
  // Get geometry index (iterative update should be done after for the dual resolution part)
  // and check grid UID
  size_t igeom;
  if (dualResolutionGridUid_ == "") {
    igeom = 0;
    if (util::getGridUid(fset) != gridUid_) {
      *infoChannel_ << "BUMP::iterativeUpdate: wrong grid UID" << std::endl;
      std::abort();
    }
  } else {
    igeom = 1;
    if (util::getGridUid(fset) != dualResolutionGridUid_) {
      *infoChannel_ << "BUMP::iterativeUpdate: wrong dual resolution grid UID" << std::endl;
      std::abort();
    }
  }

  // Print info
  *infoChannel_ << "Info     :"
                << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
  *infoChannel_ << "Info     : +++ Load member " << ie+1 << " / " << nens_[igeom] << std::endl;

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
    *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    *infoChannel_ << "Info     : +++ End of iterative update" << std::endl;
    *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
  }
}

// -----------------------------------------------------------------------------

void BUMP::writeAtlasFiles() const {
  if (bumpConf_.has("output atlas files")) {
    std::vector<eckit::LocalConfiguration> outputAtlasFilesConf
      = bumpConf_.getSubConfigurations("output atlas files");
    for (const auto & output : this->fieldsToWrite(outputAtlasFilesConf)) {
      // Get file configuration
      eckit::LocalConfiguration file = this->getFileConf(*comm_, output.first);

      // Print FieldSet norm
      *testChannel_ << "+++ Norm of output parameter " << output.second.name() << ": "
                    << util::normFieldSet(output.second, vars_, *comm_)
                    << std::endl;

      // Write FieldSet
      util::writeFieldSet(*comm_, file, output.second);
    }
  }
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> BUMP::fieldsToWrite(
  const std::vector<eckit::LocalConfiguration> confs) const {
  // Define outputs vector
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> outputs;

  // Write output fields
  for (const auto & conf : confs) {
    // Get parameter and component
    const std::string param = conf.getString("parameter");
    const int icmp = conf.getInt("component", 1);

    // Create FieldSet
    atlas::FieldSet fset;
    fset.name() = param + " - " + std::to_string(icmp);

    // Get FieldSet
    const int npar = param.size();
    const char *cpar = param.c_str();
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      bump_get_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, fset.get());
    }

    // Get file configuration
    eckit::LocalConfiguration file = this->getFileConf(*comm_, conf);

    // Add pair
    outputs.push_back(std::make_pair(file, fset));
  }

  return outputs;
}

// -----------------------------------------------------------------------------

void BUMP::runDrivers() {
  if (waitForDualResolution_) {
    waitForDualResolution_ = false;
  } else {
    for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
      *infoChannel_ << "Info     :"
                    << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    << std::endl;
      *infoChannel_ << "Info     : +++ Run drivers for BUMP instance " << (jgrid+1)
                    << " / " << keyBUMP_.size() << std::endl;
      bump_run_drivers_f90(keyBUMP_[jgrid]);
      bump_partial_dealloc_f90(keyBUMP_[jgrid]);
    }
    *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
    *infoChannel_ << "Info     : +++ End of BUMP drivers" << std::endl;
    *infoChannel_ << "Info     :"
                  << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                  << std::endl;
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbal(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbalAd(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyVbal(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyStdDev(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyStdDev(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::randomizeNicas(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_randomize_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicas(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUV(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUVAd(atlas::FieldSet & fset) const {
  for (size_t jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

size_t BUMP::getCvSize() const {
  size_t n = 0;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    n += nTmp;
  }
  return n;
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicasSqrt(const atlas::Field & cv,
                             atlas::FieldSet & fset,
                             const size_t & offset) const {
  int index = offset;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_sqrt_f90(keyBUMP_[jgrid], cv.get(), fset.get(), index);
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    index += nTmp;
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicasSqrtAd(const atlas::FieldSet & fset,
                               atlas::Field & cv,
                               const size_t & offset) const {
  int index = offset;
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_sqrt_ad_f90(keyBUMP_[jgrid], fset.get(), cv.get(), index);
    int nTmp;
    bump_get_cv_size_f90(keyBUMP_[jgrid], nTmp);
    index += nTmp;
  }
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration BUMP::getFileConf(const eckit::mpi::Comm & comm,
                                            const eckit::Configuration & conf) const {
  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm.size()));
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

  return file;
}

// -----------------------------------------------------------------------------

}  // namespace bump_lib
