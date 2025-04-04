/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Fields.h"

#ifdef ECCODES_FOUND
#include <eccodes.h>
#include <stdio.h>
#include <stdlib.h>
#endif
#include <netcdf.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "src/Geometry.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

namespace quench {

// -----------------------------------------------------------------------------

static std::vector<quench::Interpolation> interpolationsVector;

// -----------------------------------------------------------------------------

std::vector<quench::Interpolation>& Fields::interpolations() {
  return interpolationsVector;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Geometry & geom,
               const oops::Variables & vars,
               const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (auto & var : vars_) {
    // Set number of levels
    var.setLevels(geom_->levels(var.name()));

    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
    fset_.add(field);
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  // Set fields to zero
  this->zero();

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other,
               const Geometry & geom)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Check number of levels
  for (const auto & var : vars_) {
    if (geom_->levels(var.name()) != geom.levels(var.name())) {
      throw eckit::Exception("Different number of levels for variable " + var.name()
        + ", cannot interpolate", Here());
    }
  }

  if (geom_->grid() == other.geom_->grid() && geom_->halo() == other.geom_->halo()) {
    // Copy fieldset
    fset_ = util::copyFieldSet(other.fset_);
  } else {
    // Setup interpolation
    const auto & interpolation = setupGridInterpolation(*other.geom_);

    // Create fieldset
    for (const auto & var : vars_) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
      fset_.add(field);
    }

    // Set interpolation type
    for (auto field : fset_) {
      field.metadata() = other.fset_[field.name()].metadata();
      if (!field.metadata().has("interp_type")) {
        field.metadata().set("interp_type", "default");
      }
    }

    // Copy fieldset
    atlas::FieldSet fset = util::copyFieldSet(other.fset_);

    // Horizontal interpolation
    interpolation->execute(fset, fset_);
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other,
               const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (const auto & var : vars_) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
    fset_.add(field);
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata() = other.fset_[field.name()].metadata();
    if (!field.metadata().has("interp_type")) {
      field.metadata().set("interp_type", "default");
    }
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (const auto & var : vars_) {
      atlas::Field field = fset_[var.name()];
      const atlas::Field fieldOther = other.fset_[var.name()];
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        const auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = viewOther(jnode, jlevel);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

Fields::Fields(const Fields & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << classname() << "::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields and copy data
  for (const auto & var : vars_) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
    const atlas::Field fieldOther = other.fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewOther(jnode, jlevel);
        }
      }
    }
    fset_.add(field);
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata() = other.fset_[field.name()].metadata();
    if (!field.metadata().has("interp_type")) {
      field.metadata().set("interp_type", "default");
    }
  }

  oops::Log::trace() << classname() << "::Fields done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::zero() {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
    }
  }
  fset_.set_dirty(false);

  oops::Log::trace() << "Fields::zero end" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::constantValue(const double & value) {
  oops::Log::trace() << classname() << "::constantValue starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) view(jnode, jlevel) = value;
        }
      }
    }
  }
  fset_.set_dirty(false);

  oops::Log::trace() << classname() << "::constantValue end" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::constantValue(const std::vector<double> & profile) {
  oops::Log::trace() << classname() << "::constantValue starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      ASSERT(field.shape(1) == static_cast<int>(profile.size()));
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) view(jnode, jlevel) = profile[jlevel];
        }
      }
    }
  }
  fset_.set_dirty(false);

  oops::Log::trace() << classname() << "::constantValue end" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::constantValue(const eckit::Configuration & config) {
  oops::Log::trace() << "Fields::constantValue starting" << std::endl;
  for (const auto & group : config.getSubConfigurations("constant group-specific value")) {
    const std::vector<std::string> vars = group.getStringVector("variables");
    const double value = group.getDouble("constant value");
    for (const auto & var : vars_) {
      if (std::find(vars.begin(), vars.end(), var.name()) != vars.end()) {
        atlas::Field field = fset_[var.name()];
        const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
        const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
        if (field.rank() == 2) {
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(0.0);
          for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              if (gmaskView(jnode, jlevel) == 1) view(jnode, jlevel) = value;
            }
          }
        }
      }
    }
  }
  fset_.set_dirty(false);

  oops::Log::trace() << classname() << "::constantValue end" << std::endl;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    const atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewRhs(jnode, jlevel);
        }
      }
      field.set_dirty(fieldRhs.dirty());
    }
  }
  time_ = rhs.time_;

  oops::Log::trace() << classname() << "::operator= end" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator+=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  // Right-hand side fieldset
  atlas::FieldSet fsetRhs;
  if (geom_->grid() == rhs.geom_->grid() && geom_->halo() == rhs.geom_->halo()) {
    // Same geometry
    fsetRhs = util::shareFields(rhs.fset_);
  } else {
    // Interpolate
    const Fields rhsInterp(rhs, *geom_);

    // Copy fieldset
    fsetRhs = util::copyFieldSet(rhsInterp.fset_);
  }

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const atlas::Field fieldRhs = fsetRhs[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) += viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator-=(const Fields & rhs) {
  oops::Log::trace() << classname() << "::operator-= starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const atlas::Field fieldRhs = rhs.fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) -= viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Fields & Fields::operator*=(const double & zz) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) *= zz;
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void Fields::axpy(const double & zz,
                  const Fields & rhs) {
  oops::Log::trace() << classname() << "::axpy starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const atlas::Field fieldRhs = rhs.fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) += zz * viewRhs(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldRhs.dirty());
    }
  }

  oops::Log::trace() << classname() << "::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

double Fields::dot_product_with(const Fields & fld2) const {
  oops::Log::trace() << classname() << "::dot_product_with starting" << std::endl;

  double zz = 0;
  const auto ownedView = atlas::array::make_view<int, 2>(geom_->fields().field("owned"));
  for (const auto & var : vars_) {
    const atlas::Field field1 = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const atlas::Field field2 = fld2.fset_[var.name()];
    if (field1.rank() == 2) {
      const auto view1 = atlas::array::make_view<double, 2>(field1);
      const auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ownedView(jnode, 0) == 1) {
            zz += view1(jnode, jlevel)*view2(jnode, jlevel);
          }
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  oops::Log::trace() << classname() << "::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void Fields::schur_product_with(const Fields & fld2) {
  oops::Log::trace() << classname() << "::schur_product_with starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const atlas::Field field2 = fld2.fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      const auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) *= view2(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || field2.dirty());
    }
  }

  oops::Log::trace() << classname() << "::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::random() {
  oops::Log::trace() << classname() << "::random starting" << std::endl;

  fset_.clear();
  for (size_t groupIndex = 0; groupIndex < geom_->groups(); ++groupIndex) {
    // Mask and ghost points fields
    const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());

    // Total size
    size_t n = 0;
    std::vector<std::string> groupVars;
    for (const auto & var : vars_) {
      if (geom_->groupIndex(var.name()) == groupIndex) {
        groupVars.push_back(var.name());
      }
    }
    for (const auto & var : groupVars) {
      const atlas::Field field = fset_[var];
      if (field.rank() == 2) {
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) ++n;
          }
        }
      }
    }
    geom_->getComm().allReduceInPlace(n, eckit::mpi::sum());

    // Local masks
    atlas::FieldSet localMasks;
    localMasks.add(geom_->fields()[gmaskName]);
    localMasks.add(geom_->functionSpace().ghost());

    // Global masks
    atlas::FieldSet globalMasks;
    atlas::Field gmaskGlobal = geom_->functionSpace().createField<int>(
      atlas::option::name(gmaskName) | atlas::option::levels(geom_->levels(groupIndex))
      | atlas::option::global());
    globalMasks.add(gmaskGlobal);
    atlas::Field ghostGlobal = geom_->functionSpace().createField<int>(atlas::option::name("ghost")
     | atlas::option::global());
    globalMasks.add(ghostGlobal);

    // Global data
    atlas::FieldSet globalData;
    for (const auto & var : vars_) {
      if (geom_->groupIndex(var.name()) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(
          atlas::option::name(var.name())
          | atlas::option::levels(geom_->levels(var.name())) | atlas::option::global());
        globalData.add(field);
      }
    }

    // Gather masks on main processor
    if (geom_->functionSpace().type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      fs.gather(localMasks, globalMasks);
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      // NodeColumns
      if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
        fs.gather(localMasks, globalMasks);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(geom_->functionSpace());
        fs.gather(localMasks, globalMasks);
      }
    } else {
      throw eckit::NotImplemented(geom_->functionSpace().type() +
        " function space not supported yet", Here());
    }

    if (geom_->getComm().rank() == 0) {
      // Random vector
      util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);

      // Copy random values
      n = 0;
      const auto ghostView = atlas::array::make_view<int, 1>(globalMasks["ghost"]);
      for (const auto & var : vars_) {
        if (geom_->groupIndex(var.name()) == groupIndex) {
          atlas::Field field = globalData[var.name()];
          const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
          const auto gmaskView = atlas::array::make_view<int, 2>(globalMasks[gmaskName]);
          if (field.rank() == 2) {
            auto view = atlas::array::make_view<double, 2>(field);
            for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
              for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
                  view(jnode, jlevel) = rand_vec[n];
                  ++n;
                }
              }
            }
          }
        }
      }
    }

    // Local data
    atlas::FieldSet localData;
    for (const auto & var : vars_) {
      if (geom_->groupIndex(var.name()) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(
          atlas::option::name(var.name()) | atlas::option::levels(var.getLevels()));
        localData.add(field);
      }
    }

    // Scatter data from main processor
    if (geom_->functionSpace().type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      fs.scatter(globalData, localData);
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      // NodeColumns
      if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
        fs.scatter(globalData, localData);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(geom_->functionSpace());
        fs.scatter(globalData, localData);
      }
    } else {
      throw eckit::NotImplemented(geom_->functionSpace().type() +
        " function space not supported yet", Here());
    }

    // Copy data
    for (const auto & var : vars_) {
      if (geom_->groupIndex(var.name()) == groupIndex) {
        fset_.add(localData[var.name()]);
      }
    }
  }

  fset_.set_dirty();  // code is too complicated, mark dirty to be safe

  // Set duplicate points to the same value
  resetDuplicatePoints();

  oops::Log::trace() << "Fields::random done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::sqrt() {
  oops::Log::trace() << classname() << "::sqrt starting" << std::endl;
  util::sqrtFieldSet(fset_);
  oops::Log::trace() << classname() << "::sqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::dirac(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::dirac starting" << std::endl;

  if (config.has("file")) {
    // Input file
    const eckit::LocalConfiguration file(config, "file");
    this->read(file);
  } else {
    // Get dirac specifications
    std::vector<double> lon = config.getDoubleVector("lon");
    std::vector<double> lat = config.getDoubleVector("lat");
    std::vector<atlas::idx_t> level = config.getIntVector("level");
    std::vector<std::string> vars = config.getStringVector("variable");

    // Check sizes
    if (lon.size() != lat.size()) throw eckit::UserError("Inconsistent dirac specification size",
      Here());
    if (lon.size() != level.size()) throw eckit::UserError("Inconsistent dirac specification size",
      Here());
    if (lon.size() != vars.size()) throw eckit::UserError("Inconsistent dirac specification size",
      Here());

    // Build KDTree for each MPI task
    const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
    const auto ownedView = atlas::array::make_view<int, 2>(geom_->fields().field("owned"));
    const auto lonlatView = atlas::array::make_view<double, 2>(geom_->functionSpace().lonlat());
    atlas::idx_t n = 0;
    for (atlas::idx_t jnode = 0; jnode < geom_->functionSpace().size(); ++jnode) {
      if ((ghostView(jnode) == 0) && (ownedView(jnode, 0) == 1)) {
        ++n;
      }
    }
    atlas::util::IndexKDTree search;
    search.reserve(n);
    for (atlas::idx_t jnode = 0; jnode < geom_->functionSpace().size(); ++jnode) {
      if ((ghostView(jnode) == 0) && (ownedView(jnode, 0) == 1)) {
        atlas::PointLonLat pointLonLat(lonlatView(jnode, 0), lonlatView(jnode, 1));
        pointLonLat.normalise();
        atlas::PointXY point(pointLonLat);
        search.insert(point, jnode);
      }
    }
    search.build();

    // Set fields to zero
    this->zero();

    // Set dirac points
    for (size_t jdir = 0; jdir < lon.size(); ++jdir) {
      // Get field
      atlas::Field field = fset_[vars[jdir]];

      // Find MPI task
      atlas::PointLonLat pointLonLat(lon[jdir], lat[jdir]);
      pointLonLat.normalise();

      // Search nearest neighbor
      size_t index = std::numeric_limits<size_t>::max();
      double distance = std::numeric_limits<double>::max();
      bool potentialConflict = false;
      if (geom_->functionSpace().size() > 0) {
        atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 2);
        index = neighbor[0].payload();
        distance = neighbor[0].distance();
        potentialConflict = (std::abs(neighbor[0].distance()-neighbor[1].distance()) < 1.0e-12);
      }
      std::vector<double> distances(geom_->getComm().size());
      geom_->getComm().allGather(distance, distances.begin(), distances.end());
      const std::vector<double>::iterator distanceMin = std::min_element(std::begin(distances),
      std::end(distances));
      size_t sameDistanceCount = 0;
      for (size_t jj = 0; jj < geom_->getComm().size(); ++jj) {
        if (std::abs(distances[jj]-*distanceMin) < 1.0e-12) {
          ++sameDistanceCount;
        }
      }
      if (sameDistanceCount > 1) {
        throw eckit::UserError("requested dirac point exactly between two gridpoints", Here());
      }

      // Find local task
      size_t localTask(-1);
      if (geom_->getComm().rank() == 0) {
        localTask = std::distance(std::begin(distances), distanceMin);
      }
      geom_->getComm().broadcast(localTask, 0);

      if (geom_->getComm().rank() == localTask) {
        // Check potential conflict
        if (potentialConflict) {
          throw eckit::UserError("requested dirac point exactly between two gridpoints", Here());
        }

        // Add Dirac impulse
        if (field.rank() == 2) {
          auto view = atlas::array::make_view<double, 2>(field);
          view(index, level[jdir]-1) = 1.0;
        }
      }

      // Print longitude / latitude / level
      double lonDir = 0.0;
      double latDir = 0.0;
      if (geom_->getComm().rank() == localTask) {
        lonDir = lonlatView(index, 0);
        latDir = lonlatView(index, 1);
      }
      geom_->getComm().allReduceInPlace(lonDir, eckit::mpi::sum());
      geom_->getComm().allReduceInPlace(latDir, eckit::mpi::sum());
      oops::Log::info() << "Info     : Dirac point #" << jdir << ": " << lonDir << " / " << latDir
        << " / " << level[jdir] << std::endl;
    }
  }

  // Set duplicate points to the same value
  resetDuplicatePoints();

  oops::Log::trace() << classname() << "::dirac done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::diff(const Fields & x1,
                  const Fields & x2) {
  oops::Log::trace() << classname() << "::diff starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    atlas::Field fieldx1 = x1.fset_[var.name()];
    atlas::Field fieldx2 = x2.fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewx1 = atlas::array::make_view<double, 2>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 2>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) = viewx1(jnode, jlevel)-viewx2(jnode, jlevel);
          }
        }
      }
      field.set_dirty(fieldx1.dirty() || fieldx2.dirty());
    }
  }

  oops::Log::trace() << classname() << "::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::toFieldSet starting" << std::endl;

  // Share internal fieldset
  fset.clear();
  fset = util::shareFields(fset_);
  for (auto field : fset) {
    field.metadata() = fset_[field.name()].metadata();
    if (!field.metadata().has("interp_type")) {
      field.metadata().set("interp_type", "default");
    }
    field.set_dirty(fset_[field.name()].dirty());
  }

  oops::Log::trace() << classname() << "::toFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::trace() << classname() << "::fromFieldSet starting" << std::endl;

  // Check input fieldset
  ASSERT(!fset.empty());

  // Reset internal fieldset
  fset_.clear();
  fset_ = util::shareFields(fset);

  // Reset variables
  vars_ = oops::Variables(fset_.field_names());
  for (const auto & field : fset_) {
    vars_[field.name()].setLevels(field.shape(1));
  }

  // Set duplicate points to the same value
  resetDuplicatePoints();

  oops::Log::trace() << classname() << "::fromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::read(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Update variables names
  oops::Variables vars_in_file;
  for (const auto & var : vars_) {
    std::string newVar = var.name();
    for (const auto & item : geom_->alias()) {
      if (item.getString("in code") == var.name()) {
        newVar = item.getString("in file");
      }
    }
    vars_in_file.push_back({newVar, var.metaData(), var.getLevels()});
  }

  // Get IO format
  const std::string ioFormat = config.getString("format", "default");

  // Read with specified IO format
  if (ioFormat == "default") {
    // Default OOPS writer

    // Create variableSizes
    std::vector<size_t> variableSizes;
    for (const auto & var : vars_in_file) {
      variableSizes.push_back(var.getLevels());
    }

    // Update configuration
    eckit::LocalConfiguration conf(config);
    if (!conf.has("latitude south to north")) {
      conf.set("latitude south to north", geom_->latSouthToNorth());
    }

    // Read fieldset
    util::readFieldSet(geom_->getComm(),
                       geom_->functionSpace(),
                       variableSizes,
                       vars_in_file.variables(),
                       conf,
                       fset_);
  } else if (ioFormat == "grib") {
#ifdef ECCODES_FOUND
    // GRIB format

    // Build filepath
    std::string filepath = config.getString("filepath");
    if (config.has("member")) {
      std::ostringstream out;
      out << std::setfill('0') << std::setw(6) << config.getInt("member");
      filepath.append("_");
      filepath.append(out.str());
    }

    // Grib file path
    std::string gribfilepath = filepath;
    gribfilepath.append(".");
    gribfilepath.append(config.getString("grib extension", "grib2"));

    // Get levels
    std::vector<size_t> levels;
    if (!config.get("levels", levels)) {
      size_t levelMax = 0;
      for (const auto & var : vars_in_file) {
        levelMax = std::max(levelMax, geom_->levels(var.name()));
      }
      for (size_t k = 0; k < levelMax; ++k) {
        levels.push_back(k+1);
      }
    }

    // Clear local fieldset
    fset_.clear();

    // Create local fieldset
    for (const auto & var : vars_in_file) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(geom_->levels(var.name())));
      fset_.add(field);
    }

    // Initialize local fieldset
    for (auto & field : fset_) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
    }

    // Global data
    atlas::FieldSet globalData;
    for (const auto & var : vars_in_file) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(var.name())
        | atlas::option::levels(geom_->levels(var.name())) | atlas::option::global());
      globalData.add(field);
    }

    // Grib input
    if (geom_->getComm().rank() == 0) {
      oops::Log::info() << "Info     : Reading file: " << gribfilepath << std::endl;

      // Initialization
      int ret;
      codes_index* index;
      codes_handle* h;

      // Create index of file contents for cfVarName, typeOfLevel and level
      index = codes_index_new_from_file(0, gribfilepath.c_str(), "cfVarName,typeOfLevel,level",
        &ret);
      CODES_CHECK(ret, 0);

      for (const auto & var : vars_in_file) {
        // Get field view
        auto varView = atlas::array::make_view<double, 2>(globalData[var.name()]);

        // Select variable and type of level
        CODES_CHECK(codes_index_select_string(index, "cfVarName", var.name().c_str()), 0);
        CODES_CHECK(codes_index_select_string(index, "typeOfLevel", "hybrid"), 0);

        for (size_t k = 0; k < geom_->levels(var.name()); ++k) {
          // Select level
          CODES_CHECK(codes_index_select_long(index, "level", levels[k]), 0);

          // Create handle
          h = codes_handle_new_from_index(index, &ret);
          CODES_CHECK(ret, 0);

          // Print all available keys
          codes_keys_iterator *kit = codes_keys_iterator_new(h, 0, NULL);
          while (codes_keys_iterator_next(kit) == 1) {
            oops::Log::debug() << "Key: " << codes_keys_iterator_get_name(kit) << std::endl;
          }

          // Get the data size
          size_t values_len = 0;
          CODES_CHECK(codes_get_size(h, "values", &values_len), 0);

          // Allocate data
          std::vector<double> values;
          values.resize(values_len);

          // Get data
          CODES_CHECK(codes_get_double_array(h, "values", values.data(), &values_len), 0);

          // Copy data to FieldSet
          for (size_t jnode = 0; jnode < values_len; ++jnode) {
            varView(jnode, k) = values[jnode];
          }

          // Delete handle
          CODES_CHECK(codes_handle_delete(h), 0);
        }
      }

      // Check number of levels
      h = codes_handle_new_from_index(index, &ret);
      if (ret == 0) {
        throw eckit::Exception("Mismatch between level numbers in file and geometry", Here());
      }

      // Delete index
      codes_index_delete(index);
    }

    // Scatter data from main processor
    if (geom_->functionSpace().type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      // NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    }

    fset_.set_dirty();  // code is too complicated, mark dirty to be safe
#else
    throw eckit::UserError("ECCODES not available", Here());
#endif
  } else if (ioFormat == "arome") {
    // AROME data at NetCDF format (converted from epygram)

    // Build filepath
    std::string filepath = config.getString("filepath");
    if (config.has("member")) {
      std::ostringstream out;
      out << std::setfill('0') << std::setw(6) << config.getInt("member");
      filepath.append("_");
      filepath.append(out.str());
    }

    // NetCDF file path
    std::string ncFilePath = filepath + ".nc";

    // Clear local fieldset
    fset_.clear();

    // Create local fieldset
    for (size_t jvar = 0; jvar < vars_in_file.size(); ++jvar) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(vars_in_file[jvar].name()) |
        atlas::option::levels(vars_in_file[jvar].getLevels()));
      fset_.add(field);
    }

    // Initialize local fieldset
    for (auto & field : fset_) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(0.0);
    }

    // NetCDF IDs
    size_t nVarLev = 0;
    std::vector<std::string> varLevName;
    for (size_t jvar = 0; jvar < vars_in_file.size(); ++jvar) {
      for (int k = 0; k < vars_in_file[jvar].getLevels(); ++k) {
        const std::string level = std::to_string(k+40);
        varLevName.push_back("S" + std::string(3-level.length(), '0') + level
          + vars_in_file[jvar].name());
        ++nVarLev;
      }
    }
    int ncid, retval, var_id[nVarLev];

    // Global data
    atlas::FieldSet globalData;
    for (size_t jvar = 0; jvar < vars_in_file.size(); ++jvar) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(vars_in_file[jvar].name())
        | atlas::option::levels(vars_in_file[jvar].getLevels()) | atlas::option::global());
      globalData.add(field);
    }

    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();

      oops::Log::info() << "Info     : Reading file: " << ncFilePath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

      // Get variables
      for (size_t jVarLev = 0; jVarLev < nVarLev; ++jVarLev) {
        if ((retval = nc_inq_varid(ncid, varLevName[jVarLev].c_str(), &var_id[jVarLev]))) {
          ERR(retval, varLevName[jVarLev]);
        }
      }

      size_t iVarLev = 0;
      for (size_t jvar = 0; jvar < vars_in_file.size(); ++jvar) {
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_in_file[jvar].name()]);
        for (int k = 0; k < vars_in_file[jvar].getLevels(); ++k) {
          // Read data
          std::vector<double> zvar(ny * nx);
          if ((retval = nc_get_var_double(ncid, var_id[iVarLev], zvar.data()))) {
            ERR(retval, varLevName[iVarLev]);
          }
          ++iVarLev;

          // Copy data
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
              atlas::gidx_t gidx = grid.index(i, ny-1-j);
              varView(gidx, k) = zvar[j*nx + i];
            }
          }
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Scatter data from main processor
    fs.scatter(globalData, fset_);

    fset_.set_dirty();  // code is too complicated, mark dirty to be safe
  } else {
    throw eckit::UserError("Unknown I/O format", Here());
  }

  // Rename fields
  for (auto & field : fset_) {
    for (const auto & item : geom_->alias()) {
      if (item.getString("in file") == field.name()) {
        field.rename(item.getString("in code"));
      }
    }
  }

  // Set interpolation type
  for (auto field : fset_) {
    field.metadata().set("interp_type", "default");
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::write(const eckit::Configuration & config) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Copy fieldset
  atlas::FieldSet fset = util::copyFieldSet(fset_);

  // Rename fields
  for (auto & field : fset) {
    for (const auto & item : geom_->alias()) {
      if (item.getString("in code") == field.name()) {
        field.rename(item.getString("in file"));
      }
    }
  }

  // Get IO format
  const std::string ioFormat = config.getString("format", "default");

  // Write with specified IO format
  if (ioFormat == "default") {
    // Update configuration
    eckit::LocalConfiguration conf(config);
    if (!conf.has("latitude south to north")) {
      conf.set("latitude south to north", geom_->latSouthToNorth());
    }

    // Default OOPS writer
    util::writeFieldSet(geom_->getComm(), conf, fset);
  } else if (ioFormat == "grib") {
    // GRIB format
    throw eckit::NotImplemented("GRIB output not implemented yet", Here());
  } else if (ioFormat == "arome") {
    // Default OOPS writer
    util::writeFieldSet(geom_->getComm(), config, fset);
  } else {
    throw eckit::UserError("Unknown I/O format", Here());
  }

  if (geom_->mesh().generated() && config.getBool("write gmsh", false)) {
    // GMSH file path
    std::string gmshfilepath = config.getString("filepath");;
    gmshfilepath.append(".msh");
    oops::Log::info() << "Info     : Writing file: " << gmshfilepath << std::endl;

    // GMSH configuration
    const auto gmshConfig =
    atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", true) |
    atlas::util::Config("info", true);
    atlas::output::Gmsh gmsh(gmshfilepath, gmshConfig);

     // Write GMSH
    gmsh.write(geom_->mesh());
    gmsh.write(fset, fset[0].functionspace());
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

double Fields::norm() const {
  oops::Log::trace() << classname() << "::norm" << std::endl;
  return util::normFieldSet(fset_, vars_.variables(), geom_->getComm());
}

// -----------------------------------------------------------------------------

void Fields::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << std::endl;
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix << "  Geometry: " << geom_->grid().name() << " [" << geom_->grid().size() << "]"
    << std::endl;
  os << prefix << "  Fields:";
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_) {
    os << std::endl;
    double zz = 0.0;
    atlas::Field field = fset_[var.name()];
    const std::string gmaskName = "gmask_" + std::to_string(geom_->groupIndex(var.name()));
    const auto gmaskView = atlas::array::make_view<int, 2>(geom_->fields()[gmaskName]);
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zz += view(jnode, jlevel)*view(jnode, jlevel);
          }
        }
      }
    }
    geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
    zz = std::sqrt(zz);
    os << prefix << "    " << var.name() << ": " << zz;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t Fields::serialSize() const {
  oops::Log::trace() << classname() << "::serialSize starting" << std::endl;

  size_t nn = 0;
  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    if (field.rank() == 2) {
      nn += field.shape(0)*field.shape(1);
    }
  }

  oops::Log::trace() << classname() << "::serialSize done" << std::endl;
  return nn;
}

// -----------------------------------------------------------------------------

void Fields::serialize(std::vector<double> & vect)  const {
  oops::Log::trace() << classname() << "::serialize starting" << std::endl;

  for (const auto & var : vars_) {
    const atlas::Field field = fset_[var.name()];
    if (field.rank() == 2) {
      const auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          vect.push_back(view(jnode, jlevel));
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void Fields::deserialize(const std::vector<double> & vect,
                         size_t & index) {
  oops::Log::trace() << classname() << "::deserialize starting" << std::endl;

  for (const auto & var : vars_) {
    atlas::Field field = fset_[var.name()];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = vect[index];
          ++index;
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<Interpolation>::iterator Fields::setupGridInterpolation(const Geometry & srcGeom)
  const {
  oops::Log::trace() << classname() << "::setupGridInterpolation starting" << std::endl;

  // Get geometry UIDs (grid + "_" + paritioner)
  const std::string srcGeomUid = srcGeom.grid().uid() + "_" + srcGeom.partitioner().type();
  const std::string geomUid = geom_->grid().uid() + "_" + geom_->partitioner().type();

  // Compare with existing UIDs
  for (auto it = interpolations().begin(); it != interpolations().end(); ++it) {
    if ((it->srcUid() == srcGeomUid) && (it->tgtUid() == geomUid)) {
      oops::Log::trace() << classname() << "::setupGridInterpolation done" << std::endl;
      return it;
    }
  }

  // Create interpolation
  Interpolation interpolation(srcGeom,
                              srcGeomUid,
                              geom_->grid(),
                              geom_->functionSpace(),
                              geomUid);

  // Insert new interpolation
  interpolations().push_back(interpolation);

  oops::Log::trace() << classname() << "::setupGridInterpolation done" << std::endl;
  return std::prev(interpolations().end());
}

// -----------------------------------------------------------------------------

void Fields::resetDuplicatePoints() {
  oops::Log::trace() << classname() << "::resetDuplicatePoints starting" << std::endl;

  if (geom_->duplicatePoints()) {
    if (geom_->gridType() == "regular_lonlat") {
      // Deal with poles
      for (auto field_internal : fset_) {
        // Get first longitude value
        atlas::functionspace::StructuredColumns fs(field_internal.functionspace());
        atlas::StructuredGrid grid = fs.grid();
        auto view = atlas::array::make_view<double, 2>(field_internal);
        auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
        auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
        std::vector<double> north(field_internal.shape(1), 0.0);
        std::vector<double> south(field_internal.shape(1), 0.0);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (view_i(jnode) == 1) {
              if (view_j(jnode) == 1) {
                for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                  north[jlevel] = view(jnode, jlevel);
                }
              }
              if (view_j(jnode) == grid.ny()) {
                for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                  south[jlevel] = view(jnode, jlevel);
                }
              }
            }
          }
        }

        // Reduce
        geom_->getComm().allReduceInPlace(north.begin(), north.end(), eckit::mpi::sum());
        geom_->getComm().allReduceInPlace(south.begin(), south.end(), eckit::mpi::sum());

        // Copy value
        for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
          for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (view_j(jnode) == 1) {
              for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                view(jnode, jlevel) = north[jlevel];
              }
            }
            if (view_j(jnode) == grid.ny()) {
              for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
                view(jnode, jlevel) = south[jlevel];
              }
            }
          }
        }
      }
    } else {
      throw eckit::NotImplemented("duplicate points not supported for this grid", Here());
    }
  }

  oops::Log::trace() << classname() << "::resetDuplicatePoints done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
