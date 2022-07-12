/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Fields.h"

#include <netcdf.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
#include "atlas/meshgenerator/detail/CubedSphereDualMeshGenerator.h"
#include "atlas/meshgenerator/detail/CubedSphereMeshGenerator.h"
#endif
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "quench/Geometry.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(1);}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Fields::Fields(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time):
  geom_(new Geometry(geom)), vars_(vars), time_(time)
{
  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    fset_.add(field);
  }

  // Set fields to zero
  this->zero();
}
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const Geometry & geom):
  geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_)
{
  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    fset_.add(field);
  }

  // Copy - TODO(Benjamin): interpolate
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldOther = other.fset_[vars_[jvar]];
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
// -----------------------------------------------------------------------------
Fields::Fields(const Fields & other, const bool copy):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  // Reset ATLAS fieldset
  fset_ = atlas::FieldSet();

  // Create fields
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    fset_.add(field);
  }

  // Set fields to zero
  this->zero();

  // Copy if necessary
  if (copy) {
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field = fset_[vars_[jvar]];
      atlas::Field fieldOther = other.fset_[vars_[jvar]];
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
  fset_ = atlas::FieldSet();

  // Create fields and copy data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = geom_->functionSpace().createField<double>(
      atlas::option::name(vars_[jvar]) | atlas::option::levels(geom_->levels()));
    atlas::Field fieldOther = other.fset_[vars_[jvar]];
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
    fset_.add(field);
  }
}
// -----------------------------------------------------------------------------
void Fields::zero() {
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldRhs = rhs.fset_[vars_[jvar]];
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldRhs = rhs.fset_[vars_[jvar]];
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) += viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldRhs = rhs.fset_[vars_[jvar]];
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) -= viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) *= zz;
      }
    } else if (field.rank() == 2) {
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldRhs = rhs.fset_[vars_[jvar]];
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewRhs = atlas::array::make_view<double, 1>(fieldRhs);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) *= zz;
        view(jnode) += viewRhs(jnode);
      }
    } else if (field.rank() == 2) {
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
  if (geom_->functionSpace().type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field1 = fset_[vars_[jvar]];
      atlas::Field field2 = fld2.fset_[vars_[jvar]];
      if (field1.rank() == 1) {
        auto view1 = atlas::array::make_view<double, 1>(field1);
        auto view2 = atlas::array::make_view<double, 1>(field2);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            zz += view1(jnode)*view2(jnode);
          }
        }
      } else if (field1.rank() == 2) {
        auto view1 = atlas::array::make_view<double, 2>(field1);
        auto view2 = atlas::array::make_view<double, 2>(field2);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
              zz += view1(jnode, jlevel)*view2(jnode, jlevel);
            }
          }
        }
      }
    }
  } else if (geom_->functionSpace().type() == "NodeColumns") {
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field1 = fset_[vars_[jvar]];
      atlas::Field field2 = fld2.fset_[vars_[jvar]];
      if (field1.rank() == 1) {
        auto view1 = atlas::array::make_view<double, 1>(field1);
        auto view2 = atlas::array::make_view<double, 1>(field2);
        for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
          zz += view1(jnode)*view2(jnode);
        }
      } else if (field1.rank() == 2) {
        auto view1 = atlas::array::make_view<double, 2>(field1);
        auto view2 = atlas::array::make_view<double, 2>(field2);
        for (atlas::idx_t jnode = 0; jnode < field1.shape(0); ++jnode) {
          for (atlas::idx_t jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          zz += view1(jnode, jlevel)*view2(jnode, jlevel);
          }
        }
      }
    }
  }
  this->geom_->getComm().allReduceInPlace(zz, eckit::mpi::sum());
  return zz;
}
// -----------------------------------------------------------------------------
void Fields::schur_product_with(const Fields & dx) {
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldDx = dx.fset_[vars_[jvar]];
    if (field.rank() == 1) {
      auto view = atlas::array::make_view<double, 1>(field);
      auto viewDx = atlas::array::make_view<double, 1>(fieldDx);
      for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
        view(jnode) *= viewDx(jnode);
      }
    } else if (field.rank() == 2) {
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
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
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
    atlas::Field field = fset_[vars_[jvar]];
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
    if (index[jdir] <= 0 || index[jdir] > geom_->grid().size()) {
      ABORT("dirac index is too large");
    }
    if (!vars_.has(variable[jdir])) {
      ABORT("dirac variable is wrong");
    }

    if (geom_->functionSpace().type() == "StructuredColumns") {
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      atlas::Field field_gi = fs.global_index();
      auto view_gi = atlas::array::make_view<atlas::gidx_t, 1>(field_gi);
      atlas::Field field_var = fset_[variable[jdir]];
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
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      if (geom_->grid().name().substr(0, 2).compare("CS") == 0) {
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
        atlas::functionspace::CubedSphereNodeColumns fs(geom_->functionSpace());
        atlas::Field field_gi = fs.global_index();
        auto view_gi = atlas::array::make_view<atlas::gidx_t, 1>(field_gi);
        atlas::Field field_var = fset_[variable[jdir]];
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
#endif
      }
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::diff(const Fields & x1, const Fields & x2) {
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = fset_[vars_[jvar]];
    atlas::Field fieldx1 = x1.fset_[vars_[jvar]];
    atlas::Field fieldx2 = x2.fset_[vars_[jvar]];
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
void Fields::toFieldSet(atlas::FieldSet & fset) const {
  for (auto var : vars_.variables()) {
    if (fset_.has_field(var)) {
      fset->add(fset_[var]);
      atlas::Field field_input = fset_[var];
      atlas::Field field_local = fset[var];
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
      ABORT("Variable " + var + " not in source fieldset");
    }
  }
}
// -----------------------------------------------------------------------------
void Fields::fromFieldSet(const atlas::FieldSet & fset) {
  for (auto var : vars_.variables()) {
    if (fset_.has_field(var)) {
      if (fset.has_field(var)) {
        atlas::Field field_input = fset_[var];
        atlas::Field field_local = fset[var];
        if (field_input != field_local) {
          auto view_input = atlas::array::make_view<double, 2>(field_input);
          auto view_local = atlas::array::make_view<double, 2>(field_local);
          for (atlas::idx_t jnode = 0; jnode < field_input.shape(0); ++jnode) {
            for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
              view_input(jnode, jlevel) = view_local(jnode, jlevel);
            }
          }
        }

        if (geom_->gridConfig().getString("type") == "regular_lonlat") {
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
                  north[jlevel] = view(jnode,jlevel);
                }
              }
              if ((view_j(jnode) == grid.ny())  && (view_i(jnode) == 1)) {
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
                  south[jlevel] = view(jnode,jlevel);
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
                  view(jnode,jlevel) = north[jlevel];
                }
              }
              if (view_j(jnode) == grid.ny()) {
                for (atlas::idx_t jlevel = 0; jlevel < field_input.shape(1); ++jlevel) {
                  view(jnode,jlevel) = south[jlevel];
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

  // NetCDF input
  if (geom_->functionSpace().type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    // Create global data fieldset
    atlas::FieldSet globalData;
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      atlas::Field field = fs.createField<double>(atlas::option::name(vars_[jvar])
        | atlas::option::levels(geom_->levels()) | atlas::option::global());
      globalData.add(field);
    }

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get first field
      atlas::Field field = globalData.field(0);

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();
      atlas::idx_t nz = field.levels();

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
        double zvar[nx][ny][nz];
        if ((retval = nc_get_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);

        // Copy data
        atlas::Field field = globalData[vars_[jvar]];
        auto varView = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t k = 0; k < nz; ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
              atlas::gidx_t gidx = grid.index(i, j);
              varView(gidx, k) = zvar[i][j][k];
            }
          }
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }

    // Scatter data from main processor
    fs.scatter(globalData, fset_);
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

  // NetCDF output
  if (geom_->functionSpace().type() == "StructuredColumns") {
    atlas::functionspace::StructuredColumns fs(geom_->functionSpace());

    // Create local coordinates fieldset
    atlas::FieldSet localCoordinates;
    atlas::Field lonLocal = fs.createField<double>(atlas::option::name("lon"));
    localCoordinates.add(lonLocal);
    atlas::Field latLocal = fs.createField<double>(atlas::option::name("lat"));
    localCoordinates.add(latLocal);
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    auto lonView = atlas::array::make_view<double, 1>(lonLocal);
    auto latView = atlas::array::make_view<double, 1>(latLocal);
    for (atlas::idx_t jnode = 0; jnode < lonLocal.shape(0); ++jnode) {
       lonView(jnode) = lonlatView(jnode, 0);
       latView(jnode) = lonlatView(jnode, 1);
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
    fs.gather(fset_, globalData);

    if (geom_->getComm().rank() == 0) {
      // Get grid
      atlas::StructuredGrid grid = fs.grid();

      // Get first field
      atlas::Field field = globalData.field(0);

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();
      atlas::idx_t nz = field.levels();

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

      // Define variables
      d3D_id[0] = nx_id;
      d3D_id[1] = ny_id;
      d3D_id[2] = nz_id;
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        if ((retval = nc_def_var(ncid, vars_[jvar].c_str(), NC_DOUBLE, 3, d3D_id,
          &var_id[jvar]))) ERR(retval);
      }

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval);

      // Copy coordinates
      double zlon[nx][ny];
      double zlat[nx][ny];
      auto lonView = atlas::array::make_view<double, 1>(lonGlobal);
      auto latView = atlas::array::make_view<double, 1>(latGlobal);
      for (atlas::idx_t j = 0; j < ny; ++j) {
        for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
          atlas::gidx_t gidx = grid.index(i, j);
          zlon[i][j] = lonView(gidx);
          zlat[i][j] = latView(gidx);
        }
      }

      // Write coordinates
      if ((retval = nc_put_var_double(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_put_var_double(ncid, lat_id, &zlat[0][0]))) ERR(retval);

      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        // Copy data
        atlas::Field field = globalData[vars_[jvar]];
        auto varView = atlas::array::make_view<double, 2>(field);
        double zvar[nx][ny][nz];
        for (atlas::idx_t k = 0; k < nz; ++k) {
          for (atlas::idx_t j = 0; j < ny; ++j) {
            for (atlas::idx_t i = 0; i < grid.nx(j); ++i) {
              atlas::gidx_t gidx = grid.index(i, j);
              zvar[i][j][k] = varView(gidx, k);
            }
          }
        }

        // Write data
        if ((retval = nc_put_var_double(ncid, var_id[jvar], &zvar[0][0][0]))) ERR(retval);
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }
  }

  // GMSH file path
  std::string gmshfilepath = filepath;
  gmshfilepath.append(".msh");
  oops::Log::info() << "Writing file: " << gmshfilepath << std::endl;

  // GMSH configuration
  const auto gmshConfig =
  atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", true) |
  atlas::util::Config("info", true);
  atlas::output::Gmsh gmsh(gmshfilepath, gmshConfig);

  if (geom_->functionSpace().type() == "StructuredColumns") {
    const auto meshGen = atlas::MeshGenerator("structured");
    const auto mesh = atlas::Mesh(meshGen.generate(geom_->grid()));
    gmsh.write(mesh);
    gmsh.write(fset_, fset_[0].functionspace());
  } else if (geom_->functionSpace()->type() == "NodeColumns") {
    if (geom_->grid().name().substr(0, 2).compare("CS") == 0) {
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
      const auto meshConfig = atlas::util::Config("partitioner", "cubedsphere");
      const auto meshGen = atlas::MeshGenerator("cubedsphere_dual", meshConfig);
      const auto mesh = atlas::Mesh(meshGen.generate(geom_->grid()));
      gmsh.write(mesh);
      gmsh.write(fset_, fset_[0].functionspace());
#endif
    }
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
    if (geom_->functionSpace().type() == "StructuredColumns") {
      atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
      atlas::Field field = fset_[vars_[jvar]];
      if (field.rank() == 1) {
        auto view = atlas::array::make_view<double, 1>(field);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            zz += view(jnode)*view(jnode);
          }
        }
      } else if (field.rank() == 2) {
        auto view = atlas::array::make_view<double, 2>(field);
        for (atlas::idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
          for (atlas::idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              zz += view(jnode, jlevel)*view(jnode, jlevel);
            }
          }
        }
      }
    } else if (geom_->functionSpace().type() == "NodeColumns") {
      atlas::Field field = fset_[vars_[jvar]];
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
    zz = sqrt(zz);
    os << "  " << vars_[jvar] << ": " << zz << std::endl;
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
