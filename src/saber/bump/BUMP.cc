/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/BUMP.h"

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

BUMP::BUMP(const eckit::mpi::Comm & comm,
           const atlas::FunctionSpace & functionSpace1,
           const atlas::FieldSet & extraFields1,
           const std::vector<size_t> & variableSizes1,
           const oops::Variables & activeVars,
           const BUMPParameters & params,
           const std::vector<atlas::FieldSet> & fsetVec1,
           const size_t & ens1_ne_in,
           const atlas::FunctionSpace & functionSpace2,
           const atlas::FieldSet & extraFields2,
           const std::vector<atlas::FieldSet> & fsetVec2,
           const size_t & ens2_ne_in) :
  params_(params), activeVars_(activeVars), keyBUMP_(), membersConfig1_(), membersConfig2_(),
  activeVarsPerGrid_() {
  oops::Log::trace() << "BUMP::BUMP construction starting" << std::endl;

  // If testing is activated, replace _MPI_ and _OMP_ patterns
  if (params_.general.value().testing.value()) {
    // Convert to eckit configuration
    eckit::LocalConfiguration fullConfig;
    params_.serialize(fullConfig);

    // Get number of MPI tasks and OpenMP threads
    std::string mpi(std::to_string(comm.size()));
    std::string omp("1");
#ifdef _OPENMP
    # pragma omp parallel
    {
      omp = std::to_string(omp_get_num_threads());
    }
#endif
    oops::Log::info() << "Info     : MPI tasks:      " << mpi << std::endl;
    oops::Log::info() << "Info     : OpenMP threads: " << omp << std::endl;
    oops::Log::info() << fullConfig << std::endl;

    // Replace patterns
    util::seekAndReplace(fullConfig, "_MPI_", mpi);
    util::seekAndReplace(fullConfig, "_OMP_", omp);

    // Convert back to parameters
    params_.deserialize(fullConfig);
  }

  // Get ensemble 1 size
  int ens1_ne = ens1_ne_in;
  const boost::optional<eckit::LocalConfiguration> &ensembleConfig1 = params_.ensemble1.value();
  if (ensembleConfig1 != boost::none) {
    // Abort if both "members" and "members from template" are specified
    if (ensembleConfig1->has("members") && ensembleConfig1->has("members from template"))
      ABORT("BUMP: both members and members from template are specified");

    if (ensembleConfig1->has("members")) {
      // Explicit members
      ensembleConfig1->get("members", membersConfig1_);
      ens1_ne = membersConfig1_.size();
    } else if (ensembleConfig1->has("members from template")) {
      // Templated members
      eckit::LocalConfiguration templateConfig;
      ensembleConfig1->get("members from template", templateConfig);
      eckit::LocalConfiguration membersTemplate;
      templateConfig.get("template", membersTemplate);
      std::string pattern;
      templateConfig.get("pattern", pattern);
      templateConfig.get("nmembers", ens1_ne);
      int start = 1;
      if (templateConfig.has("start")) {
        templateConfig.get("start", start);
      }
      std::vector<int> except;
      if (templateConfig.has("except")) {
        templateConfig.get("except", except);
      }
      int zpad = 0;
      if (templateConfig.has("zero padding")) {
        templateConfig.get("zero padding", zpad);
      }
      int count = start;
      for (int ie=0; ie < ens1_ne; ++ie) {
        while (std::count(except.begin(), except.end(), count)) {
          count += 1;
        }
        eckit::LocalConfiguration memberConfig(membersTemplate);
        util::seekAndReplace(memberConfig, pattern, count, zpad);
        membersConfig1_.push_back(memberConfig);
        count += 1;
      }
    } else {
      ABORT("BUMP: ensemble 1 not specified");
    }
  }

  // Get ensemble 2 size
  int ens2_ne = ens2_ne_in;
  const boost::optional<eckit::LocalConfiguration> &ensembleConfig2 = params_.ensemble2.value();
  if (ensembleConfig2 != boost::none) {
    // Abort if both "members" and "members from template" are specified
    if (ensembleConfig2->has("members") && ensembleConfig2->has("members from template"))
      ABORT("BUMP: both members and members from template are specified");

    if (ensembleConfig2->has("members")) {
      // Explicit members
      ensembleConfig2->get("members", membersConfig2_);
      ens2_ne = membersConfig2_.size();
    } else if (ensembleConfig2->has("members from template")) {
      // Templated members
      eckit::LocalConfiguration templateConfig;
      ensembleConfig2->get("members from template", templateConfig);
      eckit::LocalConfiguration membersTemplate;
      templateConfig.get("template", membersTemplate);
      std::string pattern;
      templateConfig.get("pattern", pattern);
      templateConfig.get("nmembers", ens2_ne);
      int start = 1;
      if (templateConfig.has("start")) {
        templateConfig.get("start", start);
      }
      std::vector<int> except;
      if (templateConfig.has("except")) {
        templateConfig.get("except", except);
      }
      int zpad = 0;
      if (templateConfig.has("zero padding")) {
        templateConfig.get("zero padding", zpad);
      }
      int count = start;
      for (int ie=0; ie < ens2_ne; ++ie) {
        while (std::count(except.begin(), except.end(), count)) {
          count += 1;
        }
        eckit::LocalConfiguration memberConfig(membersTemplate);
        util::seekAndReplace(memberConfig, pattern, count, zpad);
        membersConfig2_.push_back(memberConfig);
        count += 1;
      }
    } else {
      ABORT("BUMP: ensemble 2 not specified");
    }
  }

  // Copy universe radius
  atlas::FieldSet universe_radius = atlas::FieldSet();
  for (const auto & fset : fsetVec1) {
    if (fset.name() == "universe radius") {
      for (const auto & field : fset) {
        universe_radius.add(field);
      }
    }
  }

  // Initialize configuration
  eckit::LocalConfiguration conf(params_.toConfiguration());

  // Update ensemble sizes
  conf.set("ensemble sizes.total ensemble size", ens1_ne);
  conf.set("ensemble sizes.total lowres ensemble size", ens2_ne);

  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get the grids configuration from input configuration and complete it
  if (conf.has("grids")) {
    // Get grids from input configuration
    conf.get("grids", grids);
    ASSERT(grids.size() > 0);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Check grids number
  ASSERT(grids.size() > 0);

  // Loop over grids
  for (auto & grid : grids) {
    // Merge conf into grid
    grid = util::mergeConfigs(grid, conf);

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    grid.get("model.variables", vars_str);
    if (vars_str.size() == 0) {
      vars_str = activeVars_.variables();
      grid.set("model.variables", vars_str);
    }

    // Save variables for each grid
    const oops::Variables gridVars(vars_str);
    activeVarsPerGrid_.push_back(gridVars);

    // Get the required number of levels add it to the grid configuration
    int nl0 = 0;
    for (size_t jvar = 0; jvar < activeVars_.size(); ++jvar) {
      std::string var = activeVars_[jvar];
      if (gridVars.has(var)) {
        int nl0_tmp = static_cast<int>(variableSizes1[jvar]);
        nl0 = std::max(nl0, std::max(nl0_tmp, 1));
      }
    }
    grid.set("model.nl0", nl0);

    // Add level index for 2D fields (first or last, first by default)
    if (!grid.has("model.lev2d")) {
      grid.set("model.lev2d", "first");
    }

    // Create BUMP instance
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, &comm, functionSpace1.get(), extraFields1.get(),
                    grid, universe_radius.get());
    keyBUMP_.push_back(keyBUMP);

    // Second geometry
    if (params_.drivers.value().compute_cov2.value()
      || params_.drivers.value().compute_cor2.value()
      || params_.drivers.value().compute_loc2.value()) {
      bump_second_geometry_f90(keyBUMP, functionSpace2.get(), extraFields2.get());
    }
  }

  // Get max number of components (geometry 1)
  bool parametersToPass1 = false;
  size_t ncmp1 = 1;
  for (const auto & fset : fsetVec1) {
    if (fset.name() != "universe radius" && fset.name() != "ensemble member") {
      parametersToPass1 = true;
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        size_t component = std::stoi(fset.name().substr(pos+2));
        ncmp1 = std::max(ncmp1, component);
      }
    }
  }
  if (parametersToPass1) {
    // Set number of components (geometry 1)
    this->setNcmp(1, ncmp1);

    // Set parameters (geometry 1)
    for (const auto & fset : fsetVec1) {
      if (fset.name() != "universe radius" && fset.name() != "ensemble member") {
        int component = 1;
        std::string name = fset.name();
        size_t pos = fset.name().find("::");
        if (pos != std::string::npos) {
          // Get component (geometry 1)
          component = std::stoi(fset.name().substr(pos+2));
          name = fset.name().substr(0, pos);
        }
        this->setParameter(name, component, fset);
      }
    }
  }

  // Get max number of components (geometry 2)
  bool parametersToPass2 = false;
  size_t ncmp2 = 1;
  for (const auto & fset : fsetVec2) {
    if (fset.name() != "universe radius" && fset.name() != "ensemble member") {
      parametersToPass2 = true;
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        size_t component = std::stoi(fset.name().substr(pos+2));
        ncmp2 = std::max(ncmp2, component);
      }
    }
  }
  if (parametersToPass2) {
    // Set number of components (geometry 2)
    this->setNcmp(2, ncmp2);

    // Set parameters (geometry 2)
    for (const auto & fset : fsetVec2) {
      if (fset.name() != "universe radius" && fset.name() != "ensemble member") {
        int component = 1;
        std::string name = fset.name();
        size_t pos = fset.name().find("::");
        if (pos != std::string::npos) {
          // Get component (geometry 2)
          component = std::stoi(fset.name().substr(pos+2));
          name = fset.name().substr(0, pos);
        }
        this->setParameter(name, component, fset);
      }
    }
  }

  oops::Log::trace() << "BUMP:BUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

BUMP::~BUMP() {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    if (keyBUMP_[jgrid] > 0) bump_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

void BUMP::addMember(const atlas::FieldSet & fset, const int & ie,
                     const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_add_member_f90(keyBUMP_[jgrid], fset.get(), ie+1, iens);
  }
}

// -----------------------------------------------------------------------------

void BUMP::updateVbalCov(const atlas::FieldSet & fset, const int & ie) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_vbal_cov_f90(keyBUMP_[jgrid], fset.get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

void BUMP::updateVar(const atlas::FieldSet & fset, const int & ie) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_var_f90(keyBUMP_[jgrid], fset.get(), ie+1);
  }
}

// -----------------------------------------------------------------------------

void BUMP::updateMom(const atlas::FieldSet & fset, const int & ie,
                            const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_update_mom_f90(keyBUMP_[jgrid], fset.get(), ie+1, iens);
  }
}

// -----------------------------------------------------------------------------

void BUMP::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_run_drivers_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbal(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyVbal(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyVbalAd(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyStdDev(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::inverseMultiplyStdDev(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::randomizeNicas(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_randomize_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyNicas(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUV(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::multiplyPsiChiToUVAd(atlas::FieldSet & fset) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_psichi_to_uv_ad_f90(keyBUMP_[jgrid], fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::getParameter(const std::string & param, const int & icmp,
  const int & igeom, atlas::FieldSet & fset) const {
  const int npar = param.size();
  const char *cpar = param.c_str();
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_get_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, igeom, fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::setNcmp(const int & igeom, const int & ncmp) const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_set_ncmp_f90(keyBUMP_[jgrid], igeom, ncmp);
  }
}

// -----------------------------------------------------------------------------

void BUMP::setParameter(const std::string & param, const int & icmp,
  const atlas::FieldSet & fset) const {
  const int npar = param.size();
  const char *cpar = param.c_str();
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_set_parameter_f90(keyBUMP_[jgrid], npar, cpar, icmp, fset.get());
  }
}

// -----------------------------------------------------------------------------

void BUMP::partialDealloc() const {
  for (unsigned int jgrid = 0; jgrid < keyBUMP_.size(); ++jgrid) {
    bump_partial_dealloc_f90(keyBUMP_[jgrid]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
