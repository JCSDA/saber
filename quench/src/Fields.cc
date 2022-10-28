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
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"

#include "src/Geometry.h"

#include "saber/src/saber/interpolation/AtlasInterpWrapper.h"

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
      atlas::option::name(var) | atlas::option::levels(geom_->variableSize(var)));
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
  if (geom_->levels() != geom.levels()) ABORT("different number of levels, cannot interpolate");

  for (const auto & var : vars_.variables()) {
    // Create field
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->variableSize(var)));
    fset_.add(field);
  }

  // Interpolate
  saber::interpolation::AtlasInterpWrapper interp(other.geom_->partitioner(),
    other.geom_->functionSpace(), geom.grid(), geom.functionSpace());
  interp.execute(other.fset_, fset_);

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
      atlas::option::name(var) | atlas::option::levels(geom_->variableSize(var)));
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
      atlas::option::name(var) | atlas::option::levels(geom_->variableSize(var)));
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
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) view(jnode, jlevel) = 0.0;
        }
      }
    }
  }
  oops::Log::trace() << "Fields::zero end" << std::endl;
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
    }
  }
  time_ = rhs.time_;
  oops::Log::trace() << "Fields::operator=(const Fields & rhs) end" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator+=(const Fields & rhs) {
  oops::Log::trace() << "Fields::operator+=(const Fields & rhs) starting" << std::endl;
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
            view(jnode, jlevel) += viewRhs(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::operator+=(const Fields & rhs) done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator-=(const Fields & rhs) {
  oops::Log::trace() << "Fields::operator-=(const Fields & rhs) starting" << std::endl;
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
            view(jnode, jlevel) -= viewRhs(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::operator-=(const Fields & rhs) done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator*=(const double & zz) {
  oops::Log::trace() << "Fields::operator*=(const Fields & rhs) starting" << std::endl;
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
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
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
            view(jnode, jlevel) *= zz;
            view(jnode, jlevel) += viewRhs(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::axpy done" << std::endl;
}
// -----------------------------------------------------------------------------
double Fields::dot_product_with(const Fields & fld2) const {
  oops::Log::trace() << "Fields::dot_product_with starting" << std::endl;
  double zz = 0;
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    atlas::Field field1 = fset_[var];
    atlas::Field field2 = fld2.fset_[var];
    if (field1.rank() == 2) {
      auto view1 = atlas::array::make_view<double, 2>(field1);
      auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1 && ghostView(jnode) == 0) {
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
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldDx = dx.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewDx = atlas::array::make_view<double, 2>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
            view(jnode, jlevel) *= viewDx(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::schur_product_with done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::random() {
  oops::Log::trace() << "Fields::random starting" << std::endl;
  // Total size
  size_t n = 0;
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1 && ghostView(jnode) == 0) ++n;
        }
      }
    }
  }
  geom_->getComm().allReduceInPlace(n, eckit::mpi::sum());

  // Local masks
  atlas::FieldSet localMasks;
  localMasks.add(geom_->extraFields().field("gmask"));
  localMasks.add(geom_->functionSpace().ghost());

  // Global masks
  atlas::FieldSet globalMasks;
  atlas::Field gmaskGlobal = geom_->functionSpace().createField<int>(atlas::option::name("gmask")
    | atlas::option::levels(geom_->levels()) | atlas::option::global());
  globalMasks.add(gmaskGlobal);
  atlas::Field ghostGlobal = geom_->functionSpace().createField<int>(atlas::option::name("ghost")
    | atlas::option::global());
  globalMasks.add(ghostGlobal);

  // Global data
  atlas::FieldSet globalData;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
      | atlas::option::levels(geom_->variableSize(var)) | atlas::option::global());
    globalData.add(field);
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
    util::NormalDistribution<double>rand_vec(n, 0.0, 1.0, 1);

    // Copy random values
    n = 0;
    auto gmaskView = atlas::array::make_view<int, 2>(globalMasks.field("gmask"));
    auto ghostView = atlas::array::make_view<int, 1>(globalMasks.field("ghost"));
    for (const auto & var : vars_.variables()) {
      atlas::Field field = globalData[var];
      if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1 && ghostView(jnode) == 0) {
              view(jnode, jlevel) = rand_vec[n];
              ++n;
            }
          }
        }
      }
    }
  }

  // Scatter data from main processor
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
    fs.scatter(globalData, fset_);
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    // NodeColumns
    if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    }
  } else {
    ABORT(geom_->functionSpace().type() + " function space not supported yet");
  }
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
  auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  auto lonlatView = atlas::array::make_view<double, 2>(geom_->functionSpace().lonlat());
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
    atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);
    size_t index(neighbor[0].payload());
    double distance(neighbor[0].distance());
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
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  for (const auto & var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldx1 = x1.fset_[var];
    atlas::Field fieldx2 = x2.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewx1 = atlas::array::make_view<double, 2>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 2>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1) {
            view(jnode, jlevel) = viewx1(jnode, jlevel)-viewx2(jnode, jlevel);
          }
        }
      }
    }
  }
  oops::Log::trace() << "Fields::diff done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::trace() << "Fields::toFieldSet starting" << std::endl;
  for (auto var : vars_.variables()) {
    if (fset_.has(var)) {
      fset->add(fset_[var]);
      atlas::Field field_input = fset_[var];
      atlas::Field field_local = fset[var];
      auto view_input = atlas::array::make_view<double, 2>(field_input);
      auto view_local = atlas::array::make_view<double, 2>(field_local);
      for (atlas::idx_t jnode = 0; jnode < field_input.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
          view_local(jnode, jlevel) = view_input(jnode, jlevel);
        }
      }
    } else {
      ABORT("Variable " + var + " not in source fieldset");
    }
  }
  oops::Log::trace() << "Fields::toFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::trace() << "Fields::fromFieldSet starting" << std::endl;
  // Get ghost points mask
  auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());

  // Copy values (excluding halo points)
  for (auto var : vars_.variables()) {
    if (fset_.has(var)) {
      if (fset.has(var)) {
        atlas::Field field_input = fset_[var];
        atlas::Field field_local = fset[var];
        auto view_input = atlas::array::make_view<double, 2>(field_input);
        auto view_local = atlas::array::make_view<double, 2>(field_local);
        for (atlas::idx_t jnode = 0; jnode < field_input.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
            if (ghostView(jnode) == 0) view_input(jnode, jlevel) = view_local(jnode, jlevel);
          }
        }

        if (geom_->gridType() == "regular_lonlat") {
          // Copy poles points
          atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
          atlas::StructuredGrid grid = fs.grid();
          auto view = atlas::array::make_view<double, 2>(field_input);
          auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
          auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
          std::vector<double> north(field_input.shape(1), 0.0);
          std::vector<double> south(field_input.shape(1), 0.0);
          for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
            for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
              atlas::idx_t jnode = fs.index(i, j);
              if ((view_j(jnode) == 1)  && (view_i(jnode) == 1)) {
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
                  north[jlevel] = view(jnode, jlevel);
                }
              }
              if ((view_j(jnode) == grid.ny())  && (view_i(jnode) == 1)) {
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
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
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
                  view(jnode, jlevel) = north[jlevel];
                }
              }
              if (view_j(jnode) == grid.ny()) {
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
                  view(jnode, jlevel) = south[jlevel];
                }
              }
            }
          }
        }
      } else {
        ABORT("Variable " + var + " not in destination fieldset");
      }
    } else {
      ABORT("Variable " + var + " not in source fieldset");
    }
  }
  oops::Log::trace() << "Fields::fromFieldSet done" << std::endl;
}
// -----------------------------------------------------------------------------
void Fields::read(const eckit::Configuration & config) {
  // Filepath
  std::string filepath = config.getString("filepath");
  if (config.has("member")) {
    std::ostringstream out;
    out << std::setfill('0') << std::setw(6) << config.getInt("member");
    filepath.append("_");
    filepath.append(out.str());
  }

  // Global data
  atlas::FieldSet globalData;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
      | atlas::option::levels(geom_->variableSize(var)) | atlas::option::global());
    globalData.add(field);
  }

  // NetCDF input
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();

      // NetCDF IDs
      int ncid, retval, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".");
      ncfilepath.append(config.getString("netcdf extension", "nc"));
      oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get variables
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_inq_varid(ncid, vars_[jvar].c_str(), &var_id[jvar]))) ERR(retval);
      }

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Read data
        double zvar[geom_->variableSize(vars_[jvar])][ny][nx];
        if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);

        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        for (size_t k = 0; k < geom_->variableSize(vars_[jvar]); ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
              atlas::gidx_t gidx = grid.index(i, ny-1-j);
              varView(gidx, k) = zvar[k][j][i];
            }
          }
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    // NodeColumns
    atlas::idx_t nb_nodes;
    if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    }

    if (geom_->getComm().rank() == 0) {
      // NetCDF IDs
      int ncid, retval, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get variables
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_inq_varid(ncid, vars_[jvar].c_str(), &var_id[jvar]))) ERR(retval);
      }

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Read data
        double zvar[nb_nodes][geom_->variableSize(vars_[jvar])];
        if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);

        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        for (size_t k = 0; k < geom_->variableSize(vars_[jvar]); ++k) {
          for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
            varView(i, k) = zvar[i][k];
          }
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  } else {
    ABORT(geom_->functionSpace().type() + " function space not supported yet");
  }

  // Scatter data from main processor
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
    fs.scatter(globalData, fset_);
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    // NodeColumns
    if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());
      fs.scatter(globalData, fset_);
    }
  } else {
    ABORT(geom_->functionSpace().type() + " function space not supported yet");
  }
}
// -----------------------------------------------------------------------------
void Fields::write(const eckit::Configuration & config) const {
  // Filepath
  std::string filepath = config.getString("filepath");
  if (config.has("member")) {
    std::ostringstream out;
    out << std::setfill('0') << std::setw(6) << config.getInt("member");
    filepath.append("_");
    filepath.append(out.str());
  }

  // Missing value
  const double msvalr = util::missingValue(double());
  const int msvali = util::missingValue(int());

  // Local coordinates
  atlas::FieldSet localCoordinates;
  atlas::Field lonLocal = geom_->functionSpace().createField<double>(atlas::option::name("lon"));
  localCoordinates.add(lonLocal);
  atlas::Field latLocal = geom_->functionSpace().createField<double>(atlas::option::name("lat"));
  localCoordinates.add(latLocal);
  localCoordinates.add(geom_->extraFields().field("gmask"));
  auto lonViewLocal = atlas::array::make_view<double, 1>(lonLocal);
  auto latViewLocal = atlas::array::make_view<double, 1>(latLocal);

  // Global coordinates
  atlas::FieldSet globalCoordinates;
  atlas::Field lonGlobal = geom_->functionSpace().createField<double>(atlas::option::name("lon")
    | atlas::option::global());
  globalCoordinates.add(lonGlobal);
  atlas::Field latGlobal = geom_->functionSpace().createField<double>(atlas::option::name("lat")
    | atlas::option::global());
  globalCoordinates.add(latGlobal);
  atlas::Field gmaskGlobal = geom_->functionSpace().createField<int>(atlas::option::name("gmask")
    | atlas::option::levels(geom_->levels()) | atlas::option::global());
  globalCoordinates.add(gmaskGlobal);
  auto lonViewGlobal = atlas::array::make_view<double, 1>(lonGlobal);
  auto latViewGlobal = atlas::array::make_view<double, 1>(latGlobal);
  auto gmaskViewGlobal = atlas::array::make_view<int, 2>(gmaskGlobal);

  // Global data
  atlas::FieldSet globalData;
  for (const auto & var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(atlas::option::name(var)
      | atlas::option::levels(geom_->variableSize(var)) | atlas::option::global());
    globalData.add(field);
  }

  // NetCDF output
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    // Local coordinates
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
       lonViewLocal(jnode) = lonlatView(jnode, 0);
       latViewLocal(jnode) = lonlatView(jnode, 1);
    }

    // Gather coordinates on main processor
    fs.gather(localCoordinates, globalCoordinates);

    // Gather data on main processor
    fs.gather(fset_, globalData);

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();
      atlas::idx_t nz = geom_->levels();

      // NetCDF IDs
      int ncid, retval, nx_id, ny_id, nz_id, d2D_id[2], d3D_id[3],
        lon_id, lat_id, gmask_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

      // Create NetCDF file
      if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

      // Create dimensions
      if ((retval = nc_def_dim(ncid, "nx", nx, &nx_id))) ERR(retval);
      if ((retval = nc_def_dim(ncid, "ny", ny, &ny_id))) ERR(retval);
      if ((retval = nc_def_dim(ncid, "nz", nz, &nz_id))) ERR(retval);
      d2D_id[0] = ny_id;
      d2D_id[1] = nx_id;
      d3D_id[0] = nz_id;
      d3D_id[1] = ny_id;
      d3D_id[2] = nx_id;

      // Define coordinates
      if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 2, d2D_id, &lon_id))) ERR(retval);
      if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 2, d2D_id, &lat_id))) ERR(retval);
      if ((retval = nc_def_var(ncid, "gmask", NC_INT, 3, d3D_id, &gmask_id))) ERR(retval);
      if ((retval = nc_put_att_double(ncid, lon_id, "_FillValue", NC_DOUBLE, 1, &msvalr)))
        ERR(retval);
      if ((retval = nc_put_att_double(ncid, lat_id, "_FillValue", NC_DOUBLE, 1, &msvalr)))
        ERR(retval);
      if ((retval = nc_put_att_int(ncid, gmask_id, "_FillValue", NC_INT, 1, &msvali)))
        ERR(retval);

      // Define variables
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if (geom_->variableSize(vars_[jvar]) == geom_->levels()) {
           if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 3, d3D_id,
          &var_id[jvar]))) ERR(retval);
        } else {
           if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
          &var_id[jvar]))) ERR(retval);
        }
        if ((retval = nc_put_att_double(ncid, var_id[jvar], "_FillValue", NC_DOUBLE, 1, &msvalr)))
          ERR(retval);
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[ny][nx];
      double zlat[ny][nx];
      int zgmask[nz][ny][nx];
      for (atlas::idx_t j = 0; j < ny; ++j) {
        for (atlas::idx_t i = 0; i < nx; ++i) {
          zlon[j][i] = msvalr;
          zlat[j][i] = msvalr;
          for (atlas::idx_t k = 0; k < nz; ++k) {
            zgmask[k][j][i] = msvali;
          }
        }
        for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
          atlas::gidx_t gidx = grid.index(i, ny-1-j);
          zlon[j][i] = lonViewGlobal(gidx);
          zlat[j][i] = latViewGlobal(gidx);
          for (atlas::idx_t k = 0; k < nz; ++k) {
            zgmask[k][j][i] = gmaskViewGlobal(gidx, k);
          }
        }
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);
      if ((retval = nc_put_var_int(ncid, gmask_id, &zgmask[0][0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        double zvar[geom_->variableSize(vars_[jvar])][ny][nx];
        for (size_t k = 0; k < geom_->variableSize(vars_[jvar]); ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < nx; ++i) {
              zvar[k][j][i] = msvalr;
            }
            for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
              atlas::gidx_t gidx = grid.index(i, ny-1-j);
              if (gmaskViewGlobal(gidx, k) == 1) {
                zvar[k][j][i] = varView(gidx, k);
              }
            }
          }
        }

        // Write data
        if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    // NodeColumns
    atlas::idx_t nb_nodes;
    if (geom_->grid().name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());

      // Local coordinates
      auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
      for (atlas::idx_t jnode = 0; jnode < fs.lonlat().shape(0); ++jnode) {
         lonViewLocal(jnode) = lonlatView(jnode, 0);
         latViewLocal(jnode) = lonlatView(jnode, 1);
      }

      // Gather coordinates on main processor
      fs.gather(localCoordinates, globalCoordinates);

      // Gather data on main processor
      fs.gather(fset_, globalData);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());

      // Local coordinates
      auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
      for (atlas::idx_t jnode = 0; jnode < fs.lonlat().shape(0); ++jnode) {
         lonViewLocal(jnode) = lonlatView(jnode, 0);
         latViewLocal(jnode) = lonlatView(jnode, 1);
      }

      // Gather coordinates on main processor
      fs.gather(localCoordinates, globalCoordinates);

      // Gather data on main processor
      fs.gather(fset_, globalData);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    }

    if (geom_->getComm().rank() == 0) {
      // Get number of levels
      atlas::idx_t nz = geom_->levels();

      // NetCDF IDs
      int ncid, retval, nb_nodes_id, nz_id, d1D_id[1], d2D_id[2],
        lon_id, lat_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

      // Create NetCDF file
      if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

      // Create dimensions
      if ((retval = nc_def_dim(ncid, "nb_nodes", nb_nodes, &nb_nodes_id))) ERR(retval);
      if ((retval = nc_def_dim(ncid, "nz", nz, &nz_id))) ERR(retval);

      // Define coordinates
      d1D_id[0] = nb_nodes_id;
      if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, d1D_id, &lon_id))) ERR(retval);
      if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, d1D_id, &lat_id))) ERR(retval);

      // Define variables
      d2D_id[0] = nb_nodes_id;
      d2D_id[1] = nz_id;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if (geom_->variableSize(vars_[jvar]) == geom_->levels()) {
          if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
            &var_id[jvar]))) ERR(retval);
        } else {
          if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 1, d1D_id,
            &var_id[jvar]))) ERR(retval);
        }
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[nb_nodes][1];
      double zlat[nb_nodes][1];
      for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
        zlon[i][0] = lonViewGlobal(i);
        zlat[i][0] = latViewGlobal(i);
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        double zvar[nb_nodes][geom_->variableSize(vars_[jvar])];
        for (size_t k = 0; k < geom_->variableSize(vars_[jvar]); ++k) {
          for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
            zvar[i][k] = varView(i, k);
          }
        }

        // Write data
        if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  } else if (geom_->functionSpace().type() == "PointCloud") {
    // PointCloud
    atlas::functionspace::PointCloud fs(geom_->functionSpace());

    if (geom_->getComm().rank() == 0) {
      // Get sizes
      atlas::idx_t nlocs = fs.size();

      // NetCDF IDs
      int ncid, retval, nlocs_id, nz_id[vars_.size()], d1D_id[1], d2D_id[2],
        lon_id, lat_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

      // Create NetCDF file
      if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

      // Create dimensions
      if ((retval = nc_def_dim(ncid, "nlocs", nlocs, &nlocs_id))) ERR(retval);

      // Define coordinates
      d1D_id[0] = nlocs_id;
      if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 1, d1D_id, &lon_id))) ERR(retval);
      if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 1, d1D_id, &lat_id))) ERR(retval);

      // Define variables
      d2D_id[0] = nlocs_id;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        std::string nz_nval = vars_[jvar];
        nz_nval.append("_nval");
        if ((retval = nc_def_dim(ncid, nz_nval.c_str(), geom_->variableSize(vars_[jvar]),
          &nz_id[jvar]))) ERR(retval);
        d2D_id[1] = nz_id[jvar];
        if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
          &var_id[jvar]))) ERR(retval);
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[nlocs][1];
      double zlat[nlocs][1];
      auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
      for (atlas::idx_t i = 0; i < nlocs; ++i) {
        zlon[i][0] = lonlatView(i, 0);
        zlat[i][0] = lonlatView(i, 1);
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        double zvar[nlocs][geom_->variableSize(vars_[jvar])];
        auto varView = atlas::array::make_view<double, 2>(fset_[vars_[jvar]]);
        for (size_t k = 0; k < geom_->variableSize(vars_[jvar]); ++k) {
          for (atlas::idx_t i = 0; i < nlocs; ++i) {
            zvar[i][k] = varView(i, k);
          }
        }

        // Write data
        if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  } else {
    ABORT(geom_->functionSpace().type() + " function space not supported yet");
  }

  if (geom_->mesh().generated()) {
    // GMSH file path
    std::string gmshfilepath = filepath;
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
  double zz = this->dot_product_with(*this);
  zz = sqrt(zz);
  return zz;
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
  auto gmaskView = atlas::array::make_view<int, 2>(geom_->extraFields().field("gmask"));
  auto ghostView = atlas::array::make_view<int, 1>(geom_->functionSpace().ghost());
  for (const auto & var : vars_.variables()) {
    os << std::endl;
    double zz = 0.0;
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (gmaskView(jnode, geom_->maskLevel(var, jlevel)) == 1 && ghostView(jnode) == 0) {
            zz += view(jnode, jlevel)*view(jnode, jlevel);
          }
        }
      }
    }
    if (geom_->functionSpace().type() != "PointCloud") {
      geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
    }
    zz = sqrt(zz);
    os << prefix << "  " << var << ": " << zz;
  }
}
// -----------------------------------------------------------------------------
size_t Fields::serialSize() const {
  size_t nn = 0;
  if (geom_->functionSpace().type() == "StructuredColumns") {
    nn = geom_->functionSpace().size();
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    nn = geom_->functionSpace().size();
  }
  return nn;
}
// -----------------------------------------------------------------------------
void Fields::serialize(std::vector<double> & vect)  const {
  ABORT("not implemented yet");
}
// -----------------------------------------------------------------------------
void Fields::deserialize(const std::vector<double> & vect, size_t & index) {
  ABORT("not implemented yet");
}
// -----------------------------------------------------------------------------
}  // namespace quench
