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
#include "oops/assimilation/Increment4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#if !ATLASIFIED
#include "oops/generic/UnstructuredGrid.h"
#endif
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
  typedef oops::Increment4D<MODEL> Increment4D_;

 public:
  OoBump(const Geometry_ &, const oops::Variables &, const std::vector<util::DateTime> &,
         const eckit::LocalConfiguration);
  explicit OoBump(OoBump &);
  ~OoBump();

  // C++ interfaces
  size_t getSize() {return keyOoBump_.size();}
  int getKey(int igrid) const {return keyOoBump_[igrid];}
  void clearKey() {keyOoBump_.clear();}

  // Fortran interfaces
  void addMember(Increment4D_ &, const int &) const;
  void addPseudoMember(Increment4D_ &, const int &) const;
  void removeMember(Increment4D_ &, const int &) const;
  void removePseudoMember(Increment4D_ &, const int &) const;
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
  void multiplyNicas(Increment4D_ &) const;
  void multiplyNicas(const Increment4D_ &, Increment4D_ &) const;
  void inverseMultiplyNicas(const Increment4D_ &, Increment4D_ &) const;
  void randomize(Increment_ &) const;
  void randomize(Increment4D_ &) const;
  void getParameter(const std::string &, Increment4D_ &) const;
  void setParameter(const std::string &, const Increment4D_ &) const;

  // Aliases for inversion with GMRESR
  void multiply(const Increment_ & dxi, Increment_ & dxo) const {multiplyNicas(dxi, dxo);}
  void multiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {multiplyNicas(dxi, dxo);}

 private:
  std::vector<int> keyOoBump_;
#if !ATLASIFIED
  std::unique_ptr<oops::UnstructuredGrid> ug_;
#endif
};

// -----------------------------------------------------------------------------
template<typename MODEL>
OoBump<MODEL>::OoBump(const Geometry_ & resol,
                      const oops::Variables & vars,
                      const std::vector<util::DateTime> & timeslots,
                      const eckit::LocalConfiguration conf) : keyOoBump_() {
  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get global prefix
  std::string prefix;
  conf.get("prefix", prefix);

#if ATLASIFIED
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
    if (grids[jgrid].has("varname")) {
      grids[jgrid].get("varname", vars_str);
    } else {
      for (unsigned int jvar = 0; jvar < vars.size(); ++jvar) {
        vars_str.push_back(vars[jvar]);
      }
      grids[jgrid].set("varname", vars_str);
    }
    grids[jgrid].set("nv", vars_str.size());

    // Add input timeslots to the grid configuration
    std::vector<std::string> timeslots_str;
    if (grids[jgrid].has("timeslot")) {
      grids[jgrid].get("timeslot", timeslots_str);
    } else {
      for (unsigned int jts = 0; jts < timeslots.size(); ++jts) {
        timeslots_str.push_back(timeslots[jts].toString());
      }
      grids[jgrid].set("timeslot", timeslots_str);
    }
    grids[jgrid].set("nts", timeslots_str.size());

    // Get the required number of levels add it to the grid configuration
    Increment_ dx(resol, vars, timeslots[0]);
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());
    int nl = 0;
    for (unsigned int jvar = 0; jvar < vars_str.size(); ++jvar) {
      std::string fieldname = vars_str[jvar] + '_' + timeslots_str[0];
      atlas::Field atlasField = atlasFieldSet->field(fieldname);
      nl = std::max(nl, atlasField.levels());
    }
    grids[jgrid].set("nl", nl);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }
  }
#else
  // Get the grids configuration from the unstructured grid configuration
  Increment4D_ dx(resol, vars, timeslots);
  int colocated = 1;
  if (conf.has("colocated")) {
    colocated = conf.getInt("colocated");
  }
  ug_.reset(new oops::UnstructuredGrid(colocated, timeslots.size()));
  dx.ug_coord(*ug_.get());
  ug_->defineGeometry();
  ug_->defineGrids(grids);

  // Modify grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    int grid_index;
    grids[jgrid].get("grid_index", grid_index);
    std::ostringstream ss;
    ss << std::setw(2) << std::setfill('0') << grid_index;
    grids[jgrid].set("prefix", prefix + "_" + ss.str());
    std::vector<std::string> vars_str;
    grids[jgrid].get("variables", vars_str);
    grids[jgrid].set("varname", vars_str);
    grids[jgrid].set("nv", vars_str.size());
    std::vector<std::string> timeslots_str;
    grids[jgrid].get("timeslots", timeslots_str);
    grids[jgrid].set("timeslot", timeslots_str);
    grids[jgrid].set("nts", timeslots_str.size());
  }
#endif

  // Check grids number
  ASSERT(grids.size() > 0);

  // Print configuration
  oops::Log::info() << "Configuration: " << conf << std::endl;

  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Print configuration for this grid
    oops::Log::info() << "Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Create OoBump instance
    int keyOoBump = 0;
#if ATLASIFIED
    bump_create_f90(keyOoBump, &resol.getComm(),
                    resol.atlasFunctionSpace()->get(),
                    resol.atlasFieldSet()->get(),
                    conf, grids[jgrid]);
#else
    bump_create_f90(keyOoBump, &resol.getComm(),
                    ug_->atlasFunctionSpace()->get(),
                    ug_->atlasFieldSet()->get(),
                    conf, grids[jgrid]);
#endif
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
#if !ATLASIFIED
  ug_ = std::move(other.ug_);
#endif
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
void OoBump<MODEL>::addMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_add_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 1);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::addPseudoMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_add_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 2);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::removeMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_remove_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 1);
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::removePseudoMember(Increment4D_ & dx, const int & ie) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_remove_member_f90(keyOoBump_[jgrid], atlasFieldSet->get(), ie+1, 2);
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
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
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_inv_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyVbalInvAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_vbal_inv_ad_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyStdDev(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_stddev_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyStdDevInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_stddev_inv_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
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
void OoBump<MODEL>::multiplyNicas(Increment4D_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::multiplyNicas(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_apply_nicas_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::inverseMultiplyNicas(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  oops::IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::randomize(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_randomize_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::randomize(Increment4D_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_randomize_f90(keyOoBump_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::getParameter(const std::string & param, Increment4D_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_get_parameter_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet.get()->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoBump<MODEL>::setParameter(const std::string & param, const Increment4D_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoBump_.size(); ++jgrid) {
    bump_set_parameter_f90(keyOoBump_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
}
// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_OOBUMP_H_
