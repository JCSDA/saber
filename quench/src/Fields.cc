/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Fields.h"

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
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"

#include "src/Geometry.h"

#include "saber/interpolation/AtlasInterpWrapper.h"

#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Fields::Fields(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time):
  geom_(new Geometry(geom)), vars_(vars), time_(time)
{
  oops::Log::trace() << "Fields::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (const auto & var : vars_.variables()) {
    // Create field

    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    fset_.add(field);
  }

  // Set fields to zero
  this->zero();

  oops::Log::trace() << "Fields::Fields done" << std::endl;
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const Geometry & geom):
  geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "Fields::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Check number of levels
  for (const auto & var : vars_.variables()) {
    if (geom_->levels(var) != geom.levels(var)) {
      ABORT("different number of levels for variable " + var + ", cannot interpolate");
    }
  }

  if (geom_->grid() == other.geom_->grid() && geom_->halo() == other.geom_->halo()) {
    // Copy fieldset
    fset_ = util::copyFieldSet(other.fset_);
  } else {
    // Create fieldset
    for (const auto & var : vars_.variables()) {
      atlas::Field field = geom_->functionSpace().createField<double>(
        atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
      fset_.add(field);
    }

    // Interpolate
    saber::interpolation::AtlasInterpWrapper interp(other.geom_->partitioner(),
      other.geom_->functionSpace(), geom.grid(), geom.functionSpace());
    interp.execute(other.fset_, fset_);
  }

  oops::Log::trace() << "Fields::Fields done" << std::endl;
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const bool copy):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "Fields::Fields starting" << std::endl;

  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    fset_.add(field);
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (const auto & var : vars_.variables()) {
      atlas::Field field = fset_[var];
      atlas::Field fieldOther = other.fset_[var];
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = viewOther(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::Fields done" << std::endl;
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  oops::Log::trace() << "Fields::Fields(const Fields & other) starting" << std::endl;
  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields and copy data
  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels(var)));
    atlas::Field fieldOther = other.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewOther(jnode, jlevel);
        }
      }
    }
    fset_.add(field);
  }
  oops::Log::trace() << "Fields::Fields(const Fields & other) done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::zero() {
  oops::Log::trace() << "Fields::zero starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
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
  oops::Log::trace() << "Fields::constantValue starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
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
  oops::Log::trace() << "Fields::constantValue end" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::constantValue(const eckit::Configuration & config) {
  oops::Log::trace() << "Fields::constantValue starting" << std::endl;
  for (const auto & group : config.getSubConfigurations("constant group-specific value")) {
    const std::vector<std::string> vars = group.getStringVector("variables");
    const double value = group.getDouble("constant value");
    for (const auto & var : vars_.variables()) {
      if (std::find(vars.begin(), vars.end(), var) != vars.end()) {
        const auto gmaskView = atlas::array::make_view<int, 2>(
          geom_->fields(geom_->groupIndex(var)).field("gmask"));
        atlas::Field field = fset_[var];
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
  oops::Log::trace() << "Fields::constantValue end" << std::endl;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator=(const Fields & rhs) {
  oops::Log::trace() << "Fields::operator=(const Fields & rhs) starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewRhs(jnode, jlevel);
        }
      }
      field.set_dirty(fieldRhs.dirty());
    }
  }
  time_ = rhs.time_;
  oops::Log::trace() << "Fields::operator=(const Fields & rhs) end" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator+=(const Fields & rhs) {
  oops::Log::trace() << "Fields::operator+=(const Fields & rhs) starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
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
  oops::Log::trace() << "Fields::operator+=(const Fields & rhs) done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator-=(const Fields & rhs) {
  oops::Log::trace() << "Fields::operator-=(const Fields & rhs) starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
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
  oops::Log::trace() << "Fields::operator-=(const Fields & rhs) done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator*=(const double & zz) {
  oops::Log::trace() << "Fields::operator*=(const Fields & rhs) starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
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
  oops::Log::trace() << "Fields::operator*=(const Fields & rhs) done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
void Fields::axpy(const double & zz, const Fields & rhs) {
  oops::Log::trace() << "Fields::axpy starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
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
  oops::Log::trace() << "Fields::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
double Fields::dot_product_with(const Fields & fld2) const {
  oops::Log::trace() << "Fields::dot_product_with starting" << std::endl;
  double zz = 0;
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field1 = fset_[var];
    atlas::Field field2 = fld2.fset_[var];
    if (field1.rank() == 2) {
      auto view1 = atlas::array::make_view<double, 2>(field1);
      auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) {
            zz += view1(jnode, jlevel)*view2(jnode, jlevel);
          }
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  oops::Log::trace() << "Fields::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::schur_product_with(const Fields & dx) {
  oops::Log::trace() << "Fields::schur_product_with starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
    atlas::Field fieldDx = dx.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewDx = atlas::array::make_view<double, 2>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) *= viewDx(jnode, jlevel);
          }
        }
      }
      field.set_dirty(field.dirty() || fieldDx.dirty());
    }
  }
  oops::Log::trace() << "Fields::schur_product_with done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::random() {
  oops::Log::trace() << "Fields::random starting" << std::endl;
  for (size_t groupIndex = 0; groupIndex < geom_->groups(); ++groupIndex) {
    // Total size
    size_t n = 0;
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(groupIndex).field("gmask"));
    const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        atlas::Field field = fset_[var];
        if (field.rank() == 2) {
          for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              if (gmaskView(jnode, jlevel) == 1 && ghostView(jnode) == 0) ++n;
            }
          }
        }
      }
    }
    geom_->getComm().allReduceInPlace(n, eckit::mpi::sum());

    // Local masks
    atlas::FieldSet localMasks;
    localMasks.add(geom_->fields(groupIndex).field("gmask"));
    localMasks.add(geom_->functionSpace().ghost());

    // Global masks
    atlas::FieldSet globalMasks;
    atlas::Field gmaskGlobal = geom_->functionSpace().createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(geom_->levels(groupIndex))
      | atlas::option::global());
    globalMasks.add(gmaskGlobal);
    atlas::Field ghostGlobal = geom_->functionSpace().createField<int>(atlas::option::name("ghost")
     | atlas::option::global());
    globalMasks.add(ghostGlobal);

    // Global data
    atlas::FieldSet globalData;
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels(var)) | atlas::option::global());
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
      ABORT(geom_->functionSpace().type() + " function space not supported yet");
    }

    if (geom_->getComm().rank() == 0) {
      // Random vector
      util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);

      // Copy random values
      n = 0;
      const auto gmaskView = atlas::array::make_view<int, 2>(globalMasks.field("gmask"));
      const auto ghostView = atlas::array::make_view<int, 1>(globalMasks.field("ghost"));
      for (const auto & var : vars_.variables()) {
        if (geom_->groupIndex(var) == groupIndex) {
          atlas::Field field = globalData[var];
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
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels(var)));
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
      ABORT(geom_->functionSpace().type() + " function space not supported yet");
    }

    // Copy data
    for (const auto & var : vars_.variables()) {
      if (geom_->groupIndex(var) == groupIndex) {
        fset_.add(localData.field(var));
      }
    }
  }

  fset_.set_dirty();  // code is too complicated, mark dirty to be safe

  oops::Log::trace() << "Fields::random done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::dirac(const eckit::Configuration & config) {
  oops::Log::trace() << "Fields::dirac starting" << std::endl;
  // Get dirac specifications
  std::vector<double> lon = config.getDoubleVector("lon");
  std::vector<double> lat = config.getDoubleVector("lat");
  std::vector<atlas::idx_t> level = config.getIntVector("level");
  std::vector<std::string> variable = config.getStringVector("variable");

  // Check sizes
  if (lon.size() != lat.size()) ABORT("Inconsistent dirac specification size");
  if (lon.size() != level.size()) ABORT("Inconsistent dirac specification size");
  if (lon.size() != variable.size()) ABORT("Inconsistent dirac specification size");

  // Build KDTree for each MPI task
  atlas::util::IndexKDTree search;
  search.reserve(geom_->functionSpace().size());
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  const auto lonlatView = atlas::array::make_view<double, 2>(geom_->functionSpace().lonlat());
  atlas::idx_t n = 0;
  for (atlas::idx_t jnode = 0; jnode < geom_->functionSpace().size(); ++jnode) {
    if (ghostView(jnode) == 0) {
      atlas::PointLonLat pointLonLat(lonlatView(jnode, 0), lonlatView(jnode, 1));
      pointLonLat.normalise();
      atlas::PointXY point(pointLonLat);
      search.insert(point, n++);
    }
  }
  search.build();

  // Set fields to zero
  this->zero();

  // Set dirac points
  for (size_t jdir = 0; jdir < lon.size(); ++jdir) {
    // Get field
    atlas::Field field = fset_[variable[jdir]];

    // Find MPI task
    atlas::PointLonLat pointLonLat(lon[jdir], lat[jdir]);
    pointLonLat.normalise();

    // Search nearest neighbor
    size_t index = std::numeric_limits<size_t>::max();
    double distance = std::numeric_limits<double>::max();
    if (geom_->functionSpace().size() > 0) {
      atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);
      index = neighbor[0].payload();
      distance = neighbor[0].distance();
    }
    std::vector<double> distances(geom_->getComm().size());
    geom_->getComm().gather(distance, distances, 0);

    // Find local task
    size_t localTask(-1);
    if (geom_->getComm().rank() == 0) {
      localTask = std::distance(std::begin(distances), std::min_element(std::begin(distances),
        std::end(distances)));
    }
    geom_->getComm().broadcast(localTask, 0);

    if (geom_->getComm().rank() == localTask) {
      // Add Dirac impulse
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        view(index, level[jdir]-1) = 1.0;
      }
    }
  }
  oops::Log::trace() << "Fields::dirac done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::diff(const Fields & x1, const Fields & x2) {
  oops::Log::trace() << "Fields::diff starting" << std::endl;
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    atlas::Field field = fset_[var];
    atlas::Field fieldx1 = x1.fset_[var];
    atlas::Field fieldx2 = x2.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewx1 = atlas::array::make_view<double, 2>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 2>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, jlevel) == 1) {
            view(jnode, jlevel) = viewx1(jnode, jlevel) - viewx2(jnode, jlevel);
          }
        }
      }
      field.set_dirty(fieldx1.dirty() || fieldx2.dirty());
    }
  }
  oops::Log::trace() << "Fields::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::trace() << "Fields::toFieldSet starting" << std::endl;
  // Share internal fieldset
  fset.clear();
  fset = util::shareFields(fset_);
  for (auto field_external : fset) {
    field_external.metadata().set("interp_type", "default");
  }
  oops::Log::trace() << "Fields::toFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::trace() << "Fields::fromFieldSet starting" << std::endl;

  // Reset internal fieldset
  fset_.clear();
  fset_ = util::shareFields(fset);

  if (geom_->gridType() == "regular_lonlat") {
    // Reset poles points
    for (auto field_internal : fset_) {
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
          if ((view_j(jnode) == 1)  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              north[jlevel] = view(jnode, jlevel);
            }
          }
          if ((view_j(jnode) == grid.ny())  && (view_i(jnode) == 1)) {
            for (atlas::idx_t jlevel = 0; jlevel < field_internal.shape(1); ++jlevel) {
              south[jlevel] = view(jnode, jlevel);
            }
          }
        }
      }
      geom_->getComm().allReduceInPlace(north.begin(), north.end(), eckit::mpi::sum());
      geom_->getComm().allReduceInPlace(south.begin(), south.end(), eckit::mpi::sum());
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
  }
  oops::Log::trace() << "Fields::fromFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::read(const eckit::Configuration & config) {
  // Create variableSizes
  std::vector<size_t> variableSizes;
  for (const auto & var : vars_.variables()) {
    variableSizes.push_back(geom_->levels(var));
  }

  // Copy configuration
  eckit::LocalConfiguration conf(config);

  // Read fieldset
  util::readFieldSet(geom_->getComm(),
                     geom_->functionSpace(),
                     variableSizes,
                     vars_.variables(),
                     conf,
                     fset_);

  fset_.set_dirty();
}
// -----------------------------------------------------------------------------
void Fields::write(const eckit::Configuration & config) const {
  // Write fieldset
  util::writeFieldSet(geom_->getComm(), config, fset_);

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
    gmsh.write(fset_, fset_[0].functionspace());
  }
}
// -----------------------------------------------------------------------------
double Fields::norm() const {
  return util::normFieldSet(fset_, vars_.variables(), geom_->getComm());
}
// -----------------------------------------------------------------------------
void Fields::print(std::ostream & os) const {
  os << std::endl;
  os << *geom_;
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix << "Fields:";
  const auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    const auto gmaskView = atlas::array::make_view<int, 2>(
      geom_->fields(geom_->groupIndex(var)).field("gmask"));
    os << std::endl;
    double zz = 0.0;
    atlas::Field field = fset_[var];
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
    zz = sqrt(zz);
    os << prefix << "  " << var << ": " << zz;
  }
}
// -----------------------------------------------------------------------------
size_t Fields::serialSize() const {
  size_t nn = 0;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      nn += field.shape(0)*field.shape(1);
    }
  }
  return nn;
}
// -----------------------------------------------------------------------------
void Fields::serialize(std::vector<double> & vect)  const {
  for (const auto & var : vars_.variables()) {
    const atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      const auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          vect.push_back(view(jnode, jlevel));
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::deserialize(const std::vector<double> & vect, size_t & index) {
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
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
}
// -----------------------------------------------------------------------------
}  // namespace quench
