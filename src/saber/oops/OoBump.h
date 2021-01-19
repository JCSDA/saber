/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_OOBUMP_H_
#define SABER_OOPS_OOBUMP_H_

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "saber/bump/type_bump.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
  class UnstructuredGrid;
}

namespace saber {

// -----------------------------------------------------------------------------
/// OoBump C++ interface

template<typename MODEL> class OoBump {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment<MODEL>   Increment_;

 public:
  OoBump(const Geometry_ &, const oops::Variables &, const util::DateTime &,
         const eckit::LocalConfiguration);
  explicit OoBump(OoBump &);
  ~OoBump();

  // C++ interfaces
  size_t getSize() {return keyOoBump_.size();}
  int getKey(int igrid) const {return keyOoBump_[igrid];}
  void clearKey() {keyOoBump_.clear();}

  // Fortran interfaces
  void addMember(const atlas::FieldSet & atlasFieldSet, const int &, const int &) const;
  void runDrivers() const;
  void multiplyVbal(const Increment_ &, Increment_ &) const;
  void multiplyVbalInv(const Increment_ &, Increment_ &) const;
  void multiplyVbalAd(const Increment_ &, Increment_ &) const;
  void multiplyVbalInvAd(const Increment_ &, Increment_ &) const;
  void multiplyStdDev(const Increment_ &, Increment_ &) const;
  void multiplyStdDevInv(const Increment_ &, Increment_ &) const;
  void multiplyNicas(Increment_ &) const;
  void multiplyNicas(const Increment_ &, Increment_ &) const;
  void inverseMultiplyNicas(const Increment_ &, Increment_ &) const;
  void randomize(Increment_ &) const;
  void getParameter(const std::string &, Increment_ &) const;
  void setParameter(const std::string &, const Increment_ &) const;
  void partialDealloc() const;

  // Aliases for inversion with GMRESR
  void multiply(const Increment_ & dxi, Increment_ & dxo) const {multiplyNicas(dxi, dxo);}

 private:
  std::vector<int> keyOoBump_;
};

// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::OoBump(const Geometry_ & resol,
                      const oops::Variables & vars,
                      const util::DateTime & time,
                      const eckit::LocalConfiguration conf) : keyOoBump_() {
  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get global prefix
  std::string prefix;
  conf.get("prefix", prefix);

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

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add prefix
    if (!grids[jgrid].has("prefix")) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << jgrid;
      grids[jgrid].set("prefix", prefix + "_" + ss.str());
    }

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      for (unsigned int jvar = 0; jvar < vars.size(); ++jvar) {
        vars_str.push_back(vars[jvar]);
      }
      grids[jgrid].set("variables", vars_str);
    }
    grids[jgrid].set("nv", vars_str.size());

    // Get the required number of levels add it to the grid configuration
    Increment_ dx(resol, vars, time);
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());
    int nl = 0;
    for (unsigned int jvar = 0; jvar < vars_str.size(); ++jvar) {
      atlas::Field atlasField = atlasFieldSet->field(vars_str[jvar]);
      nl = std::max(nl, atlasField.levels());
    }
    grids[jgrid].set("nl", nl);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }
  }

  // Check grids number
  ASSERT(grids.size() > 0);

  // Print configuration
  oops::Log::info() << "Configuration: " << conf << std::endl;

  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Print configuration for this grid
    oops::Log::info() << "Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Create OoBump instance
    int keyOoBump = 0;
    bump_create_f90(keyOoBump, &resol.getComm(),
                    resol.atlasFunctionSpace()->get(),
                    resol.atlasFieldSet()->get(),
                    conf, grids[jgrid]);
    keyOoBump_.push_back(keyOoBump);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::OoBump(OoBump & other) : keyOoBump_() {
  for (unsigned int jgrid = 0; jgrid < other.getSize(); ++jgrid) {
    keyOoBump_.push_back(other.getKey(jgrid));
  }
  other.clearKey();
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::~OoBump() {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    if (keyOoBump_[jgrid] > 0) bump_dealloc_f90(keyOoBump_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::addMember(const atlas::FieldSet & atlasFieldSet, const int & ie,
                               const int & iens) const {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_add_member_f90(keyOoBump_[jgrid], atlasFieldSet.get(), ie+1, iens);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_run_drivers_f90(keyOoBump_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbal(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInvAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_inv_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyStdDev(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyStdDevInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dx.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dxi.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dxo.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::inverseMultiplyNicas(const Increment_ & dxi, Increment_ & dxo) const {
  oops::IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::randomize(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_randomize_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
  dx.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::getParameter(const std::string & param, Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.setAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_get_parameter_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
  dx.fromAtlas(atlasFieldSet.get());
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::setParameter(const std::string & param, const Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
  dx.toAtlas(atlasFieldSet.get());
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_set_parameter_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::partialDealloc() const {
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_partial_dealloc_f90(keyOoBump_[jgrid]);
  }
}
// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_OOBUMP_H_
