/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Geometry.h"

#include <math.h>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/projection.h"
#include "atlas/util/Point.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Geometry::Geometry(const Parameters_ & params,
                   const eckit::mpi::Comm & comm) : comm_(comm), levels_(1) {
  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Halo
  const boost::optional<size_t> &halo = params.halo.value();
  if (halo != boost::none) {
    halo_ = *halo;
  } else {
    halo_ = 0;
  }

  // Setup grid
  gridConfig_ = params.grid.value();
  oops::Log::info() << "Grid config: " << gridConfig_ << std::endl;
  grid_ = atlas::Grid(gridConfig_);

  if (params.functionSpace.value() == "StructuredColumns") {
    // StructuredColumns function space

    // Setup partitioner
    const atlas::grid::Partitioner partitioner(params.partitioner.value());

    // Setup distribution
    const atlas::grid::Distribution distribution(grid_, partitioner);

    // Setup function space
    functionSpace_ = atlas::functionspace::StructuredColumns(grid_, distribution,
                     atlas::option::halo(halo_));
  } else if (params.functionSpace.value() == "NodeColumns") {
    if (comm_.size() == 1) {
      if (gridConfig_.getString("name").substr(0, 2).compare("CS") == 0) {
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
        mesh_ = atlas::MeshGenerator("cubedsphere_dual").generate(grid_);
        functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(mesh_);
#endif
      } else {
        mesh_ = atlas::MeshGenerator("cubedsphere_dual").generate(grid_);
        functionSpace_ = atlas::functionspace::NodeColumns(mesh_);
      }
    } else {
      ABORT(params.functionSpace.value() + " function space on multiple PEs not implemented yet");
    }
  } else {
    ABORT(params.functionSpace.value() + " function space not implemented yet");
  }

  // Number of levels
  levels_ = params.levels.value();

  // Vertical unit
  const boost::optional<std::vector<double>> &vunitParams = params.vunit.value();
  for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
    if (vunitParams != boost::none) {
      vunit_.push_back((*vunitParams)[jlevel]);
    } else {
      vunit_.push_back(static_cast<double>(jlevel+1));
    }
  }

  // Fill extra geometry fields
  extraFields_ = atlas::FieldSet();

  // Vertical unit
  atlas::Field vunit = functionSpace_.createField<double>(
    atlas::option::name("vunit") | atlas::option::levels(levels_));
  auto view = atlas::array::make_view<double, 2>(vunit);
  for (atlas::idx_t jnode = 0; jnode < vunit.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < vunit.shape(1); ++jlevel) {
       view(jnode, jlevel) = vunit_[jlevel];
    }
  }
  extraFields_->add(vunit);

  // Halo mask
  if (gridConfig_.getString("type") == "regular_lonlat") {
    atlas::functionspace::StructuredColumns fs(functionSpace_);
    atlas::StructuredGrid grid = fs.grid();
    atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
      | atlas::option::levels(0));
    auto view = atlas::array::make_view<int, 1>(hmask);
    auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
    auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
    for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
      for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
        atlas::idx_t jnode = fs.index(i, j);
        if (((view_j(jnode) == 1) || (view_j(jnode) == grid.ny())) && (view_i(jnode) != 1)) {
          view(jnode) = 0;
        } else {
          view(jnode) = 1;
        }
      }
    }

    // Add field
    extraFields_->add(hmask);
  }

  // Print summary
  this->print(oops::Log::info());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_), levels_(other.levels_),
  vunit_(other.vunit_), halo_(other.halo_) {
  // Copy grid
  gridConfig_ = other.gridConfig_;
  grid_ = atlas::Grid(gridConfig_);

  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
  // Grid name
    if (gridConfig_.getString("name").substr(0, 2).compare("CS") == 0) {
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
#endif
    } else {
      functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    }
  } else {
    ABORT(other.functionSpace_.type() + " function space not implemented yet");
  }

  // Copy extra fields
  extraFields_ = atlas::FieldSet();
  atlas::Field vunit = (*other.extraFields())["vunit"];
  extraFields_->add(vunit);
}
// -------------------------------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes(vars.size(), levels_);
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  os << "Quench geometry grid:" << std::endl;
  os << "- name: " << grid_.name() << std::endl;
  os << "- size: " << grid_.size() << std::endl;
  os << "Function space:" << std::endl;
  os << "- type: " << functionSpace_.type() << std::endl;
  os << "- size: " << functionSpace_.size() << std::endl;
  os << "- nb_partitions: " << functionSpace_.nb_partitions() << std::endl;
  os << "- halo: " << halo_ << std::endl;
  os << "Vertical levels: " << std::endl;
  os << "- number: " << levels_ << std::endl;
  os << "- vunit: " << vunit_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace quench
