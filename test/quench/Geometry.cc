/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Geometry.h"

#include <math.h>
#include <netcdf.h>

#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/projection.h"
#include "atlas/util/Point.h"
#include "atlas/util/PolygonLocator.h"
#include "atlas/util/PolygonXY.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); ABORT("NetCDF error");}

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
  const boost::optional<eckit::LocalConfiguration> &gridParams = params.grid.value();
  const boost::optional<std::string> &gridInputParams = params.gridInput.value();
  if (gridParams != boost::none) {
    oops::Log::info() << "Grid config: " << *gridParams << std::endl;
    grid_ = atlas::Grid(*gridParams);
  } else if (gridInputParams != boost::none) {
    oops::Log::info() << "Grid input file: " << *gridInputParams << std::endl;

    // Grid size
    size_t nlocs = 0;

    // NetCDF IDs
    int ncid, retval, nlocs_id, grp_id, lon_id, lat_id;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(gridInputParams->c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get nlocs
      if ((retval = nc_inq_dimid(ncid, "nlocs", &nlocs_id))) ERR(retval);
      if ((retval = nc_inq_dimlen(ncid, nlocs_id, &nlocs))) ERR(retval);
    }

    // Broadcast size
    comm_.broadcast(nlocs, 0);

    // Coordinates
    std::vector<double> xy(2*nlocs);

    if (comm_.rank() == 0) {
      // Get lon/lat
      if ((retval = nc_inq_ncid(ncid, "MetaData", &grp_id))) ERR(retval);
      if ((retval = nc_inq_varid(grp_id, "longitude", &lon_id))) ERR(retval);
      if ((retval = nc_inq_varid(grp_id, "latitude", &lat_id))) ERR(retval);

      // Read data
      double zlon[nlocs][1];
      double zlat[nlocs][1];
      if ((retval = nc_get_var_double(grp_id, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_get_var_double(grp_id, lat_id, &zlat[0][0]))) ERR(retval);

      // Copy data
      for (size_t i = 0; i < nlocs; ++i) {
        xy[2*i] = zlon[i][0];
        xy[2*i+1] = zlat[i][0];
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }

    // Broadcast coordinates
    comm_.broadcast(xy.begin(), xy.end(), 0);

    // Define grid configuration
    eckit::LocalConfiguration gridConfig;
    gridConfig.set("type", "unstructured");
    gridConfig.set("xy", xy);
    grid_ = atlas::Grid(gridConfig);
  } else {
    ABORT("Grid or grid input file required");
  }

  // Setup partitioner
  partitioner_ = atlas::grid::Partitioner(params.partitioner.value());

  // Setup distribution
  const atlas::grid::Distribution distribution(grid_, partitioner_);

  if (params.functionSpace.value() == "StructuredColumns") {
    // StructuredColumns

    // Setup function space
    functionSpace_ = atlas::functionspace::StructuredColumns(grid_, distribution,
                     atlas::option::halo(halo_));

    // Setup mesh
    mesh_ = atlas::MeshGenerator("structured").generate(grid_, partitioner_);
  } else if (params.functionSpace.value() == "NodeColumns") {
    if (comm_.size() == 1) {
      // NodeColumns
      if (grid_.name().substr(0, 2).compare("CS") == 0) {
        // CubedSphere
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
        mesh_ = atlas::MeshGenerator("cubedsphere_dual").generate(grid_);
        functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(mesh_);
#else
        ABORT("TRANS required");
#endif
      } else {
        // NodeColumns
        mesh_ = atlas::MeshGenerator("delaunay").generate(grid_);
        functionSpace_ = atlas::functionspace::NodeColumns(mesh_);
      }
    } else {
      ABORT("NodeColumns function space on multiple PEs not supported yet");
    }
  } else if (params.functionSpace.value() == "PointCloud") {
    // Setup function space
    functionSpace_ = atlas::functionspace::PointCloud(grid_);

    // Setup mesh
//    mesh_ = atlas::MeshGenerator("delaunay").generate(grid_);
  } else {
    ABORT(params.functionSpace.value() + " function space not supported yet");
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
  if (grid_.name().substr(0, 1).compare("L") == 0) {
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
  // Copy grid TODO (in header ?)
  grid_ = other.grid_;
  partitioner_ = other.partitioner_;
  mesh_ = other.mesh_;

  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().substr(0, 2).compare("CS") == 0) {
      // CubedSphere
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
#endif
    } else {
      // Other NodeColumns
      functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    }
  } else if (other.functionSpace_.type() == "PointCloud") {
      functionSpace_ = atlas::functionspace::PointCloud(other.functionSpace_);
  } else {
    ABORT(other.functionSpace_.type() + " function space not supported yet");
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
  os << "Partitioner:" << std::endl;
  os << "- type: " << partitioner_.type() << std::endl;
  os << "Function space:" << std::endl;
  os << "- type: " << functionSpace_.type() << std::endl;
  os << "- halo: " << halo_ << std::endl;
  os << "Vertical levels: " << std::endl;
  os << "- number: " << levels_ << std::endl;
  os << "- vunit: " << vunit_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace quench
