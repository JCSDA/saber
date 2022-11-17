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

  // Parameters
  const GeneralParameters general = params_.general.value().get_value_or(GeneralParameters());
  const DriversParameters drivers = params_.drivers.value().get_value_or(DriversParameters());


  // If testing is activated, replace _MPI_ and _OMP_ patterns
  const bool testing = general.testing.value().get_value_or(false);
  if (testing) {
    // Convert to eckit configuration
    eckit::LocalConfiguration fullConfig;
    params_.serialize(fullConfig);

    // Get number of MPI tasks and OpenMP threads
    std::string mpi(std::to_string(comm.size()));
    std::string omp("1");
    # pragma omp parallel
    {
      omp = std::to_string(omp_get_num_threads());
    }
    oops::Log::info() << "Info     : MPI tasks:      " << mpi << std::endl;
    oops::Log::info() << "Info     : OpenMP threads: " << omp << std::endl;

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
  atlas::FieldSet universe_rad = atlas::FieldSet();
  for (const auto & fset : fsetVec1) {
    if (fset.name() == "universe radius") {
      for (const auto & field : fset) {
        universe_rad.add(field);
      }
    }
  }

  // Initialize configuration
  eckit::LocalConfiguration conf(params_.toConfiguration());

  // Add missing value (real)
  conf.set("msvalr", util::missingValue(double()));

  // Add ensemble sizes
  if (!conf.has("ens1_ne")) conf.set("ens1_ne", ens1_ne);
  if (!conf.has("ens2_ne")) conf.set("ens2_ne", ens2_ne);

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

  // Print configuration for this grid
  oops::Log::info() << "Info     : General configuration: " << conf << std::endl;

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      vars_str = activeVars_.variables();
      grids[jgrid].set("variables", vars_str);
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
    grids[jgrid].set("nl0", nl0);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }

    // Print configuration for this grid
    oops::Log::info() << "Info     : Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Create BUMP instance
    oops::Log::info() << "Info     : Create BUMP instance " << jgrid << std::endl;
    int keyBUMP = 0;
    bump_create_f90(keyBUMP, &comm, functionSpace1.get(), extraFields1.get(),
                    conf, grids[jgrid], universe_rad.get());
    keyBUMP_.push_back(keyBUMP);

    // Second geometry
    std::string method = drivers.method.value().get_value_or("");
    if (method == "hyb-ens" || method == "hyb-rnd") {
      bump_second_geometry_f90(keyBUMP, functionSpace2.get(), extraFields2.get());
    }
  }

  // Get max number of components
  size_t ncmp1 = 1;
  for (const auto & fset : fsetVec1) {
    if (fset.name() != "universe radius") {
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        size_t component = std::stoi(fset.name().substr(pos+2));
        ncmp1 = std::max(ncmp1, component);
      }
    }
  }
  size_t ncmp2 = 1;
  for (const auto & fset : fsetVec2) {
    if (fset.name() != "universe radius") {
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        size_t component = std::stoi(fset.name().substr(pos+2));
        ncmp2 = std::max(ncmp2, component);
      }
    }
  }

  // Set parameters
  this->setNcmp(1, ncmp1);
  for (const auto & fset : fsetVec1) {
    if (fset.name() != "universe radius") {
      int component = 1;
      std::string name = fset.name();
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        // Get component
        component = std::stoi(fset.name().substr(pos+2));
        name = fset.name().substr(0, pos);
      }
      this->setParameter(name, component, fset);
    }
  }
  this->setNcmp(2, ncmp2);
  for (const auto & fset : fsetVec2) {
    if (fset.name() != "universe radius") {
      int component = 1;
      std::string name = fset.name();
      size_t pos = fset.name().find("::");
      if (pos != std::string::npos) {
        // Get component
        component = std::stoi(fset.name().substr(pos+2));
        name = fset.name().substr(0, pos);
      }
      this->setParameter(name, component, fset);
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
