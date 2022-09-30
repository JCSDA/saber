/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Fields.h"

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
#include "oops/util/Random.h"

#include "quench/Geometry.h"

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

  // Create fields
  for (const auto var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels()));
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

  // Create fields
  for (const auto var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels()));
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

  // Create fields
  for (const auto var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels()));
    fset_.add(field);
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (const auto var : vars_.variables()) {
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
  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields and copy data
  for (const auto var : vars_.variables()) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(var) | atlas::option::levels(geom_->levels()));
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
}
// -----------------------------------------------------------------------------
void Fields::zero() {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = 0.0;
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
Fields & Fields::operator=(const Fields & rhs) {
  for (const auto var : vars_.variables()) {
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
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator+=(const Fields & rhs) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) += viewRhs(jnode, jlevel);
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator-=(const Fields & rhs) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) -= viewRhs(jnode, jlevel);
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator*=(const double & zz) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= zz;
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
void Fields::axpy(const double & zz, const Fields & rhs) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldRhs = rhs.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= zz;
          view(jnode, jlevel) += viewRhs(jnode, jlevel);
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
double Fields::dot_product_with(const Fields & fld2) const {
  double zz = 0;
  atlas::Field ghost = geom_->functionSpace().ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);
  for (const auto var : vars_.variables()) {
    atlas::Field field1 = fset_[var];
    atlas::Field field2 = fld2.fset_[var];
    if (field1.rank() == 2) {
      auto view1 = atlas::array::make_view<double, 2>(field1);
      auto view2 = atlas::array::make_view<double, 2>(field2);
      for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          if (ghostView(jnode) == 0) zz += view1(jnode, jlevel)*view2(jnode, jlevel);
        }
      }
    }
  }
  this->geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::schur_product_with(const Fields & dx) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldDx = dx.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewDx = atlas::array::make_view<double, 2>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= viewDx(jnode, jlevel);
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::random() {
  // Total size
  size_t n = 0;
  atlas::Field ghost = geom_->functionSpace().ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) n += field.shape(1);
      }
    }
  }

  // Random vector
  util::NormalDistribution<double> rand_vec(n, 0.0, 1.0, 1);

  // Copy random values
  n = 0;
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) = rand_vec[n];
            ++n;
          }
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::dirac(const eckit::Configuration & config) {
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
  atlas::idx_t n{0};
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
}
// -----------------------------------------------------------------------------
void Fields::diff(const Fields & x1, const Fields & x2) {
  for (const auto var : vars_.variables()) {
    atlas::Field field = fset_[var];
    atlas::Field fieldx1 = x1.fset_[var];
    atlas::Field fieldx2 = x2.fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewx1 = atlas::array::make_view<double, 2>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 2>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewx1(jnode, jlevel)-viewx2(jnode, jlevel);
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::toFieldSet(atlas::FieldSet & fset) const {
  for (auto var : vars_.variables()) {
    if (fset_.has_field(var)) {
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
}
// -----------------------------------------------------------------------------
void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  atlas::Field ghost = geom_->functionSpace().ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);
  for (auto var : vars_.variables()) {
    if (fset_.has_field(var)) {
      if (fset.has_field(var)) {
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

  // Common object
  atlas::FieldSet globalData;

  // NetCDF input
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    // Create global data fieldset
    for (const auto var : vars_.variables()) {
      atlas::Field field = fs.createField<double>(atlas::option::name(var)
        | atlas::option::levels(geom_->levels()) | atlas::option::global());
      globalData.add(field);
    }

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();
      atlas::idx_t nz = globalData.field(0).levels();

      // NetCDF IDs
      int ncid, retval, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".");
      ncfilepath.append(config.getString("netcdf extension", "nc"));
      oops::Log::info() << "Reading file: " << ncfilepath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get variables
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_inq_varid(ncid, vars_[jvar].c_str(), &var_id[jvar]))) ERR(retval);
      }

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Read data
        double zvar[nz][ny][nx];
        if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);

        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        for (atlas::idx_t k = 0; k < nz; ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
              atlas::gidx_t gidx = grid.index(i, j);
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

      // Create global data fieldset
      for (const auto var : vars_.variables()) {
        atlas::Field field = fs.createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels()) | atlas::option::global());
        globalData.add(field);
      }

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());

      // Create global data fieldset
      for (const auto var : vars_.variables()) {
        atlas::Field field = fs.createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels()) | atlas::option::global());
        globalData.add(field);
      }

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    }

    if (geom_->getComm().rank() == 0) {
      // Get number of levels
      atlas::idx_t nz = globalData.field(0).levels();

      // NetCDF IDs
      int ncid, retval, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Reading file: " << ncfilepath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get variables
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_inq_varid(ncid, vars_[jvar].c_str(), &var_id[jvar]))) ERR(retval);
      }

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Read data
        double zvar[nb_nodes][nz];
        if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0]))) ERR(retval);

        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        for (atlas::idx_t k = 0; k < nz; ++k) {
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
  double msv(-999.0);  // TODO(Benjamin) should be missing values

  // Common objects
  atlas::FieldSet localCoordinates;
  atlas::Field lonLocal;
  atlas::Field latLocal;
  atlas::FieldSet globalCoordinates;
  atlas::Field lonGlobal;
  atlas::Field latGlobal;
  atlas::FieldSet globalData;

  // NetCDF output
  if (geom_->functionSpace().type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    // Create local coordinates fieldset
    lonLocal = fs.createField<double>(atlas::option::name("lon"));
    localCoordinates.add(lonLocal);
    latLocal = fs.createField<double>(atlas::option::name("lat"));
    localCoordinates.add(latLocal);
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    auto lonView = atlas::array::make_view<double, 1>(lonLocal);
    auto latView = atlas::array::make_view<double, 1>(latLocal);
    for (atlas::idx_t jnode = 0; jnode < lonLocal.shape(0); ++jnode) {
       lonView(jnode) = lonlatView(jnode, 0);
       latView(jnode) = lonlatView(jnode, 1);
    }

    // Create global coordinates fieldset
    lonGlobal = fs.createField<double>(atlas::option::name("lon")
      | atlas::option::global());
    globalCoordinates.add(lonGlobal);
    latGlobal = fs.createField<double>(atlas::option::name("lat")
      | atlas::option::global());
    globalCoordinates.add(latGlobal);

    // Gather coordinates on main processor
    fs.gather(localCoordinates, globalCoordinates);

    // Create global data fieldset
    for (const auto var : vars_.variables()) {
      atlas::Field field = fs.createField<double>(atlas::option::name(var)
        | atlas::option::levels(geom_->levels()) | atlas::option::global());
      globalData.add(field);
    }

    // Gather data on main processor
    fs.gather(fset_, globalData);

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();
      atlas::idx_t nz = globalData.field(0).levels();

      // NetCDF IDs
      int ncid, retval, nx_id, ny_id, nz_id, d2D_id[2], d3D_id[3],
        lon_id, lat_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

      // Create NetCDF file
      if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER, &ncid))) ERR(retval);

      // Create dimensions
      if ((retval = nc_def_dim(ncid, "nx", nx, &nx_id))) ERR(retval);
      if ((retval = nc_def_dim(ncid, "ny", ny, &ny_id))) ERR(retval);
      if ((retval = nc_def_dim(ncid, "nz", nz, &nz_id))) ERR(retval);

      // Define coordinates
      d2D_id[0] = nx_id;
      d2D_id[1] = ny_id;
      if ((retval = nc_def_var(ncid, "lon", NC_DOUBLE, 2, d2D_id, &lon_id))) ERR(retval);
      if ((retval = nc_def_var(ncid, "lat", NC_DOUBLE, 2, d2D_id, &lat_id))) ERR(retval);
      if ((retval = nc_put_att_double(ncid, lon_id, "_FillValue", NC_DOUBLE, 1, &msv)))
        ERR(retval);
      if ((retval = nc_put_att_double(ncid, lat_id, "_FillValue", NC_DOUBLE, 1, &msv)))
        ERR(retval);

      // Define variables
      d3D_id[0] = nz_id;
      d3D_id[1] = ny_id;
      d3D_id[2] = nx_id;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 3, d3D_id,
          &var_id[jvar]))) ERR(retval);
        if ((retval = nc_put_att_double(ncid, var_id[jvar], "_FillValue", NC_DOUBLE, 1, &msv)))
          ERR(retval);
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[ny][nx];
      double zlat[ny][nx];
      auto lonView = atlas::array::make_view<double, 1>(lonGlobal);
      auto latView = atlas::array::make_view<double, 1>(latGlobal);
      for (atlas::idx_t j = 0; j < ny; ++j) {
        for (atlas::idx_t i = 0; i < nx; ++i) {
          zlon[j][i] = msv;
          zlat[j][i] = msv;
        }
        for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
          atlas::gidx_t gidx = grid.index(i, j);
          zlon[j][i] = lonView(gidx);
          zlat[j][i] = latView(gidx);
        }
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        double zvar[nz][ny][nx];
        for (atlas::idx_t k = 0; k < nz; ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < nx; ++i) {
              zvar[k][j][i] = msv;
            }
            for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
              atlas::gidx_t gidx = grid.index(i, j);
              zvar[k][j][i] = varView(gidx, k);
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

      // Create local coordinates fieldset
      lonLocal = fs.createField<double>(atlas::option::name("lon"));
      localCoordinates.add(lonLocal);
      latLocal = fs.createField<double>(atlas::option::name("lat"));
      localCoordinates.add(latLocal);
      auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
      auto lonView = atlas::array::make_view<double, 1>(lonLocal);
      auto latView = atlas::array::make_view<double, 1>(latLocal);
      for (atlas::idx_t jnode = 0; jnode < lonLocal.shape(0); ++jnode) {
         lonView(jnode) = lonlatView(jnode, 0);
         latView(jnode) = lonlatView(jnode, 1);
      }

      // Create global coordinates fieldset
      lonGlobal = fs.createField<double>(atlas::option::name("lon")
        | atlas::option::global());
      globalCoordinates.add(lonGlobal);
      latGlobal = fs.createField<double>(atlas::option::name("lat")
        | atlas::option::global());
      globalCoordinates.add(latGlobal);

      // Gather coordinates on main processor
      fs.gather(localCoordinates, globalCoordinates);

      // Create global data fieldset
      for (const auto var : vars_.variables()) {
        atlas::Field field = fs.createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels()) | atlas::option::global());
        globalData.add(field);
      }

      // Gather data on main processor
      fs.gather(fset_, globalData);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    } else {
      // Other NodeColumns
      atlas::functionspace::NodeColumns fs(geom_->functionSpace());

      // Create local coordinates fieldset
      lonLocal = fs.createField<double>(atlas::option::name("lon"));
      localCoordinates.add(lonLocal);
      latLocal = fs.createField<double>(atlas::option::name("lat"));
      localCoordinates.add(latLocal);
      auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
      auto lonView = atlas::array::make_view<double, 1>(lonLocal);
      auto latView = atlas::array::make_view<double, 1>(latLocal);
      for (atlas::idx_t jnode = 0; jnode < lonLocal.shape(0); ++jnode) {
         lonView(jnode) = lonlatView(jnode, 0);
         latView(jnode) = lonlatView(jnode, 1);
      }

      // Create global coordinates fieldset
      lonGlobal = fs.createField<double>(atlas::option::name("lon")
        | atlas::option::global());
      globalCoordinates.add(lonGlobal);
      latGlobal = fs.createField<double>(atlas::option::name("lat")
        | atlas::option::global());
      globalCoordinates.add(latGlobal);

      // Gather coordinates on main processor
      fs.gather(localCoordinates, globalCoordinates);

      // Create global data fieldset
      for (const auto var : vars_.variables()) {
        atlas::Field field = fs.createField<double>(atlas::option::name(var)
          | atlas::option::levels(geom_->levels()) | atlas::option::global());
        globalData.add(field);
      }

      // Gather data on main processor
      fs.gather(fset_, globalData);

      // Get global number of nodes
      nb_nodes = fs.nb_nodes_global();
    }

    if (geom_->getComm().rank() == 0) {
      // Get number of levels
      atlas::idx_t nz = globalData.field(0).levels();

      // NetCDF IDs
      int ncid, retval, nb_nodes_id, nz_id, d1D_id[1], d2D_id[2],
        lon_id, lat_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

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
        if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 2, d2D_id,
          &var_id[jvar]))) ERR(retval);
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[nb_nodes][1];
      double zlat[nb_nodes][1];
      auto lonView = atlas::array::make_view<double, 1>(lonGlobal);
      auto latView = atlas::array::make_view<double, 1>(latGlobal);
      for (atlas::idx_t i = 0; i < nb_nodes; ++i) {
        zlon[i][0] = lonView(i);
        zlat[i][0] = latView(i);
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        auto varView = atlas::array::make_view<double, 2>(globalData[vars_[jvar]]);
        double zvar[nb_nodes][nz];
        for (atlas::idx_t k = 0; k < nz; ++k) {
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
      atlas::idx_t nz = fset_.field(0).levels();

      // NetCDF IDs
      int ncid, retval, nlocs_id, nz_id[vars_.size()], d1D_id[1], d2D_id[2],
        lon_id, lat_id, var_id[vars_.size()];

      // NetCDF file path
      std::string ncfilepath = filepath;
      ncfilepath.append(".nc");
      oops::Log::info() << "Writing file: " << ncfilepath << std::endl;

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
        if ((retval = nc_def_dim(ncid, nz_nval.c_str(), nz, &nz_id[jvar]))) ERR(retval);
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
        double zvar[nlocs][nz];
        auto varView = atlas::array::make_view<double, 2>(fset_[vars_[jvar]]);
        for (atlas::idx_t k = 0; k < nz; ++k) {
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
    oops::Log::info() << "Writing file: " << gmshfilepath << std::endl;

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
  os << "Fields:" << std::endl;

  atlas::Field ghost = geom_->functionSpace().ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);
  for (const auto var : vars_.variables()) {
    double zz = 0.0;
    atlas::Field field = fset_[var];
    if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          if (ghostView(jnode) == 0) zz += view(jnode, jlevel)*view(jnode, jlevel);
        }
      }
    }
    this->geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
    zz = sqrt(zz);
    os << "  " << var << ": " << zz << std::endl;
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
