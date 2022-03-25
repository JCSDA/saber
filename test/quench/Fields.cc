/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Fields.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "quench/Geometry.h"
#include "quench/Interface.h"

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Fields::Fields(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time):
  geom_(new Geometry(geom)), vars_(vars), time_(time)
{
  // Reset ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());

  // Create fields
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->atlasFunctionSpace()->createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    atlasFieldSet_->add(field);
  }

  // Set fields to zero
  this->zero();
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  // Reset ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());

  // Create fields
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->atlasFunctionSpace()->createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    atlasFieldSet_->add(field);
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
      atlas::Field fieldOther = other.atlasFieldSet_->field(vars_[jvar]);
      if (field.rank() == 1) {
        auto view = atlas::array::make_view<double, 1>(field);
        auto viewOther = atlas::array::make_view<double, 1>(fieldOther);
        for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
          view(jnode) = viewOther(jnode);
        }
      } else if (field.rank() == 2) {
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
}

// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  // Reset ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());

  // Create fields and copy data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->atlasFunctionSpace()->createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    atlas::Field fieldOther = other.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewOther = atlas::array::make_view<double, 1>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = viewOther(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewOther = atlas::array::make_view<double, 2>(fieldOther);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = viewOther(jnode, jlevel);
        }
      }
    }
    atlasFieldSet_->add(field);
  }
}
// -----------------------------------------------------------------------------
void Fields::zero() {
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = 0.0;
      }
    } else if (field.rank() == 2) {
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
  // Get fields and copy data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldRhs = rhs.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
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
  // Get fields and add data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldRhs = rhs.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = view(jnode)+viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = view(jnode, jlevel)+viewRhs(jnode, jlevel);
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator-=(const Fields & rhs) {
  // Get fields and subtract data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldRhs = rhs.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = view(jnode)-viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = view(jnode, jlevel)-viewRhs(jnode, jlevel);
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
Fields & Fields::operator*=(const double & zz) {
  // Get fields and add data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = view(jnode)*zz;
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = view(jnode, jlevel)*zz;
        }
      }
    }
  }
  return *this;
}
// -----------------------------------------------------------------------------
void Fields::axpy(const double & zz, const Fields & rhs) {
  // Get fields and add data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldRhs = rhs.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = zz*view(jnode)+viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewRhs = atlas::array::make_view<double, 2>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = zz*view(jnode, jlevel)+viewRhs(jnode, jlevel);
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
double Fields::dot_product_with(const Fields & fld2) const {
  double zz = 0;
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        zz += view(jnode)*view(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        zz += view(jnode, jlevel)*view(jnode, jlevel);
        }
      }
    }
  }
  this->geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::schur_product_with(const Fields & dx) {
  // Get fields and add data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldDx = dx.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewDx = atlas::array::make_view<double, 1>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = view(jnode)*viewDx(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      auto viewDx = atlas::array::make_view<double, 2>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = view(jnode, jlevel)*viewDx(jnode, jlevel);
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::random() {
  // Total size
  size_t n = 0;
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      n += field.shape(0);
    } else if (field.rank() == 2) {
      n += field.shape(0)*field.shape(1);
    }
  }

  // Random vector
  util::NormalDistribution<double>rand_vec(n, 0.0, 1.0, 1);

  // Copy random values
  n = 0;
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = rand_vec[n];
        n += 1;
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) = rand_vec[n];
          n += 1;
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::dirac(const eckit::Configuration & config) {
  // Get dirac specifications
  std::vector<atlas::gidx_t> index = config.getLongVector("index");
  std::vector<atlas::idx_t> level = config.getIntVector("level");
  std::vector<std::string> variable = config.getStringVector("variable");

  // Set fields to zero
  this->zero();

  // Set dirac points
  for (size_t jdir = 0; jdir < index.size(); ++jdir) {
    if (index[jdir] <= 0 || index[jdir] > geom_->atlasGrid()->size()) {
      ABORT("dirac index is too large");
    }
    if (!vars_.has(variable[jdir])) {
      ABORT("dirac variable is wrong");
    }

    if (geom_->atlasFunctionSpace()->type() == "StructuredColumns") {
      atlas::functionspace::StructuredColumns fs(*(geom_->atlasFunctionSpace()));
      atlas::Field field_gi = fs.global_index();
      auto view_gi = atlas::array::make_view<atlas::gidx_t, 1>(field_gi);
      atlas::Field field_var = atlasFieldSet_->field(variable[jdir]);
      if (field_var.rank() == 1) {
        auto view_var = atlas::array::make_view<double, 1>(field_var);
        for (atlas::idx_t jnode = 0; jnode < field_var.shape(0); ++jnode) {
          if (index[jdir] == view_gi(jnode)) {
            view_var(jnode) = 1.0;
          }
        }
      } else if (field_var.rank() == 2) {
        if (level[jdir] <= 0 || level[jdir] > field_var.shape(1)) {
          ABORT("dirac level is too large");
        }
        auto view_var = atlas::array::make_view<double, 2>(field_var);
        for (atlas::idx_t jnode = 0; jnode < field_var.shape(0); ++jnode) {
          if (index[jdir] == view_gi(jnode)) {
            view_var(jnode, level[jdir]-1) = 1.0;
          }
        }
      }
    } else {
      ABORT(geom_->atlasFunctionSpace()->type() + " function space not implemented yet");
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::diff(const Fields & x1, const Fields & x2) {
  // Get fields and subtract data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldx1 = x1.atlasFieldSet_->field(vars_[jvar]);
    atlas::Field fieldx2 = x2.atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewx1 = atlas::array::make_view<double, 1>(fieldx1);
      auto viewx2 = atlas::array::make_view<double, 1>(fieldx2);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) = viewx1(jnode)-viewx2(jnode);
      }
    } else if (field.rank() == 2) {
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
void Fields::setAtlas(atlas::FieldSet * afieldset) const {
  for (auto var : vars_.variables()) {
    if (atlasFieldSet_->has_field(var)) {
      afieldset->add(atlasFieldSet_->field(var));
    } else {
      ABORT("Variable " + var + " not in increment");
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::toAtlas(atlas::FieldSet * afieldset) const {
  for (auto var : vars_.variables()) {
    if (atlasFieldSet_->has_field(var)) {
      if (afieldset->has_field(var)) {
        atlas::Field field_input = atlasFieldSet_->field(var);
        atlas::Field field_local = afieldset->field(var);
        if (field_input != field_local) {
          auto view_input = atlas::array::make_view<double, 2>(field_input);
          auto view_local = atlas::array::make_view<double, 2>(field_local);
          for (atlas::idx_t jnode = 0; jnode < field_input.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
              view_local(jnode, jlevel) = view_input(jnode, jlevel);
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
void Fields::fromAtlas(atlas::FieldSet * afieldset) {
}
// -----------------------------------------------------------------------------
void Fields::read(const eckit::Configuration & config) {
  ABORT("not implemented yet");
}
// -----------------------------------------------------------------------------
void Fields::write(const eckit::Configuration & config) const {
  // Get file path
  const std::string filepath = config.getString("filepath");

  if (geom_->atlasFunctionSpace()->type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(*(geom_->atlasFunctionSpace()));

    // Create local coordinates fieldset
    atlas::FieldSet localCoordinates;
    atlas::Field lonLocal = fs.createField<double>(atlas::option::name("lon"));
    localCoordinates.add(lonLocal);
    atlas::Field latLocal = fs.createField<double>(atlas::option::name("lat"));
    localCoordinates.add(latLocal);
    auto xyView = atlas::array::make_view<double, 2>(fs.xy());
    auto lonView = atlas::array::make_view<double, 1>(lonLocal);
    auto latView = atlas::array::make_view<double, 1>(latLocal);
    for (atlas::idx_t jnode = 0; jnode < lonLocal.shape(0); ++jnode) {
       lonView(jnode) = xyView(jnode, 0);
       latView(jnode) = xyView(jnode, 1);
    }

    // Create global coordinates fieldset
    atlas::FieldSet globalCoordinates;
    atlas::Field lonGlobal = fs.createField<double>(atlas::option::name("lon")
      | atlas::option::global());
    globalCoordinates.add(lonGlobal);
    atlas::Field latGlobal = fs.createField<double>(atlas::option::name("lat")
      | atlas::option::global());
    globalCoordinates.add(latGlobal);

    // Gather coordinates on main processor
    fs.gather(localCoordinates, globalCoordinates);

    // Create global data fieldset
    atlas::FieldSet globalData;
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field = fs.createField<double>(atlas::option::name(vars_[jvar])
        | atlas::option::levels(geom_->levels()) | atlas::option::global());
      globalData.add(field);
    }

    // Gather data on main processor
    fs.gather(*atlasFieldSet_, globalData);

    if (geom_->getComm().rank() == 0) {
      if (netcdfOutput_) {
        // Call Fortran for NetCDF output
        fields_write_structuredcolumns_f90(config, vars_, globalCoordinates.get(), globalData.get());
      } else {
        // Get file path
        const std::string filepath = config.getString("filepath");

        // Write longitudes
        atlas::Field lonGlobal = globalCoordinates.field("lon");
        auto lonView = atlas::array::make_view<double, 1>(lonGlobal);
        std::string filepath_lon = filepath + "_lon";
        std::ofstream outfile_lon(filepath_lon.c_str());
        if (outfile_lon.is_open()) {
          lonView.dump(outfile_lon);
          outfile_lon.close();
        } else {
          ABORT("Fields::write: cannot open file for longitudes");
        }

        // Write latitudes
        atlas::Field latGlobal = globalCoordinates.field("lat");
        auto latView = atlas::array::make_view<double, 1>(latGlobal);
        std::string filepath_lat = filepath + "_lat";
        std::ofstream outfile_lat(filepath_lat.c_str());
        if (outfile_lat.is_open()) {
          latView.dump(outfile_lat);
          outfile_lat.close();
        } else {
          ABORT("Fields::write: cannot open file for latitudes");
        }

        for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
          // Write variable
          atlas::Field field = globalData.field(vars_[jvar]);
          auto varView = atlas::array::make_view<double, 2>(field);
          std::string filepath_var = filepath + "_" + vars_[jvar];
          std::ofstream outfile_var(filepath_var.c_str());
          if (outfile_var.is_open()) {
            varView.dump(outfile_var);
            outfile_var.close();
          } else {
            ABORT("Fields::write: cannot open file for variable" + vars_[jvar]);
          }
        }
      }
    }
  } else {
    ABORT(geom_->atlasFunctionSpace()->type() + " function space not implemented yet");
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    double zz = 0.0;
    atlas::Field field = atlasFieldSet_->field(vars_[jvar]);
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        zz += view(jnode)*view(jnode);
      }
    } else if (field.rank() == 2) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        zz += view(jnode, jlevel)*view(jnode, jlevel);
        }
      }
    }
    zz = sqrt(zz);
    os << "  " << vars_[jvar] << ": " << zz << std::endl;
  }
}
// -----------------------------------------------------------------------------
size_t Fields::serialSize() const {
  size_t nn = 0;
  if (geom_->atlasFunctionSpace()->type() == "StructuredColumns") {
    nn = geom_->atlasFunctionSpace()->size();
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
