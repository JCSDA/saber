/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Geometry.h"

#include <netcdf.h>

#include <cmath>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------

Geometry::Geometry(const Parameters_ & params,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), levels_(params.levels.value()), lev2d_(params.lev2d.value()) {
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
  const boost::optional<std::string> &iodaFile = params.iodaFile.value();
  unstructuredGrid_ = false;
  if (gridParams != boost::none) {
    oops::Log::info() << "Info     : Grid config: " << *gridParams << std::endl;
    eckit::LocalConfiguration gridParams_(*gridParams);
    if (gridParams_.has("type")) {
      std::string type = gridParams_.getString("type");
      if (type == "unstructured") {
        // Split unstructured grid among processors
        std::vector<double> xyFull = gridParams_.getDoubleVector("xy");
        size_t gridSize = xyFull.size()/2;
        size_t rank = 0;
        std::vector<double> xy;
        for (size_t jnode = 0; jnode < gridSize; ++jnode) {
          // Copy coordinates on a given task
          if (comm_.rank() == rank) {
            xy.push_back(xyFull[2*jnode]);
            xy.push_back(xyFull[2*jnode+1]);
          }

          // Update task index
          ++rank;
          if (rank == comm_.size()) rank = 0;
        }

        // Reset coordinates
        gridParams_.set("xy", xy);

        // Set flag
        unstructuredGrid_ = true;
      }
    }
    grid_ = atlas::Grid(gridParams_);
  } else if (iodaFile != boost::none) {
    oops::Log::info() << "Info     : Grid input file: " << *iodaFile << std::endl;

    // Grid size
    size_t nlocs = 0;

    // NetCDF IDs
    int ncid, retval, nlocs_id, grp_id, lon_id, lat_id;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(iodaFile->c_str(), NC_NOWRITE, &ncid))) ERR(retval);

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

  if (!unstructuredGrid_) {
    // Setup partitioner
    partitioner_ = atlas::grid::Partitioner(params.partitioner.value());

    // Setup distribution
    distribution_ = atlas::grid::Distribution(grid_, partitioner_);
  }

  if (params.functionSpace.value() == "StructuredColumns") {
    // StructuredColumns

    // Setup function space
    functionSpace_ = atlas::functionspace::StructuredColumns(grid_, distribution_,
                     atlas::option::halo(halo_));

    // Setup mesh
    mesh_ = atlas::MeshGenerator("structured").generate(grid_, partitioner_);
  } else if (params.functionSpace.value() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      mesh_ = atlas::MeshGenerator("cubedsphere_dual").generate(grid_);
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(mesh_);
    } else {
      if (comm_.size() == 1) {
        // NodeColumns
        mesh_ = atlas::MeshGenerator("delaunay").generate(grid_);
        functionSpace_ = atlas::functionspace::NodeColumns(mesh_);
      } else {
        ABORT("NodeColumns function space on multiple PEs not supported yet");
      }
    }
  } else if (params.functionSpace.value() == "PointCloud") {
    // Setup function space
    functionSpace_ = atlas::functionspace::PointCloud(grid_);
  } else {
    ABORT(params.functionSpace.value() + " function space not supported yet");
  }

  // Ghost points
  atlas::Field ghost = functionSpace_.ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);

  // Vertical unit
  const boost::optional<std::vector<double>> &vunitParams = params.vunit.value();
  for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
    if (vunitParams != boost::none) {
      vunit_.push_back((*vunitParams)[jlevel]);
    } else {
      vunit_.push_back(static_cast<double>(jlevel+1));
    }
  }

  // Default mask, set to 1 (true)
  gmask_ = functionSpace_.createField<int>(
    atlas::option::name("gmask") | atlas::option::levels(levels_));
  auto maskView = atlas::array::make_view<int, 2>(gmask_);
  for (atlas::idx_t jnode = 0; jnode < gmask_.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < gmask_.shape(1); ++jlevel) {
       maskView(jnode, jlevel) = 1;
    }
  }

  // Specific mask
  if (params.mask_type.value() == "none") {
    // No mask
  } else if (params.mask_type.value() == "sea") {
    // Lon/lat sizes
    size_t nlon = 0;
    size_t nlat = 0;

    // NetCDF IDs
    int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(params.mask_path.value().c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      // Get lon/lat sizes
      if ((retval = nc_inq_dimid(ncid, "lon", &nlon_id))) ERR(retval);
      if ((retval = nc_inq_dimid(ncid, "lat", &nlat_id))) ERR(retval);
      if ((retval = nc_inq_dimlen(ncid, nlon_id, &nlon))) ERR(retval);
      if ((retval = nc_inq_dimlen(ncid, nlat_id, &nlat))) ERR(retval);
    }

    // Broadcast lon/lat sizes
    comm_.broadcast(nlon, 0);
    comm_.broadcast(nlat, 0);

    // Coordinates and land-sea mask
    std::vector<double> lon(nlon);
    std::vector<double> lat(nlat);
    std::vector<int> lsm(nlat*nlon);

    if (comm_.rank() == 0) {
      // Get lon/lat
      if ((retval = nc_inq_varid(ncid, "lon", &lon_id))) ERR(retval);
      if ((retval = nc_inq_varid(ncid, "lat", &lat_id))) ERR(retval);
      if ((retval = nc_inq_varid(ncid, "LSMASK", &lsm_id))) ERR(retval);

      // Read data
      float zlon[nlon][1];
      float zlat[nlat][1];
      uint8_t zlsm[nlat][nlon];
      if ((retval = nc_get_var_float(ncid, lon_id, &zlon[0][0]))) ERR(retval);
      if ((retval = nc_get_var_float(ncid, lat_id, &zlat[0][0]))) ERR(retval);
      if ((retval = nc_get_var_ubyte(ncid, lsm_id, &zlsm[0][0]))) ERR(retval);

      // Copy data
      for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lon[ilon] = zlon[ilon][0];
      }
      for (size_t ilat = 0; ilat < nlat; ++ilat) {
        lat[ilat] = zlat[ilat][0];
      }
      for (size_t ilat = 0; ilat < nlat; ++ilat) {
       for (size_t ilon = 0; ilon < nlon; ++ilon) {
          lsm[ilat*nlon+ilon] = static_cast<int>(zlsm[ilat][ilon]);
        }
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }

    // Broadcast coordinates and land-sea mask
    comm_.broadcast(lon.begin(), lon.end(), 0);
    comm_.broadcast(lat.begin(), lat.end(), 0);
    comm_.broadcast(lsm.begin(), lsm.end(), 0);

    // Build KD-tree
    atlas::Geometry geometry(atlas::util::Earth::radius());
    atlas::util::IndexKDTree2D search(geometry);
    search.reserve(nlat*nlon);
    std::vector<double> lon2d;
    std::vector<double> lat2d;
    std::vector<size_t> payload2d;
    int jnode = 0;
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
      for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lon2d.push_back(lon[ilon]);
        lat2d.push_back(lat[ilat]);
        payload2d.push_back(jnode);
        ++jnode;
      }
    }
    search.build(lon2d, lat2d, payload2d);

    if (functionSpace_.type() == "StructuredColumns") {
      // StructuredColumns
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
      for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          // Find nearest neighbor
          size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
            lonlatView(jnode, 1)}).payload();

          // Ocean points for all levels
          for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
            if (lsm[nn] == 0) {
               maskView(jnode, jlevel) = 1;
             } else {
               maskView(jnode, jlevel) = 0;
             }
           }

          // Ocean + small islands for the 2D level
          if (lsm[nn] == 3) {
            if (lev2d_ == "first") {
              maskView(jnode, 0) = 1;
            } else if (lev2d_ == "last") {
              maskView(jnode, levels_-1) = 1;
            } else {
              ABORT("wrong lev2d value (first or last)");
            }
          }
        }
      }
    } else {
      ABORT(params.mask_type.value() + " mask not supported for " + functionSpace_.type() + " yet");
    }
  } else {
    ABORT("Wrong mask type");
  }

  // Mask size
  gmaskSize_ = 0.0;
  double domainSize = 0.0;
  for (atlas::idx_t jnode = 0; jnode < gmask_.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < gmask_.shape(1); ++jlevel) {
      if (ghostView(jnode) == 0) {
        if (maskView(jnode, jlevel) == 1) {
          gmaskSize_ += 1.0;
        }
        domainSize += 1.0;
      }
    }
  }
  comm_.allReduceInPlace(gmaskSize_, eckit::mpi::sum());
  comm_.allReduceInPlace(domainSize, eckit::mpi::sum());
  if (domainSize > 0.0) {
    gmaskSize_ = gmaskSize_/domainSize;
  }

  // Fill extra geometry fields
  extraFields_ = atlas::FieldSet();

  // Vertical unit
  atlas::Field vunit = functionSpace_.createField<double>(
    atlas::option::name("vunit") | atlas::option::levels(levels_));
  auto vunitView = atlas::array::make_view<double, 2>(vunit);
  for (atlas::idx_t jnode = 0; jnode < vunit.shape(0); ++jnode) {
    for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
       vunitView(jnode, jlevel) = vunit_[jlevel];
    }
  }
  extraFields_->add(vunit);

  // Geographical mask
  extraFields_->add(gmask_);

  // Halo mask
  if (grid_.name().compare(0, 1, std::string{"L"}) == 0) {
    // Regular lonlat grid
    atlas::functionspace::StructuredColumns fs(functionSpace_);
    atlas::StructuredGrid grid = fs.grid();
    atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
      | atlas::option::levels(1));
    auto hmaskView = atlas::array::make_view<int, 2>(hmask);
    auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
    auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
    for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
      for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
        atlas::idx_t jnode = fs.index(i, j);
        if (((view_j(jnode) == 1) || (view_j(jnode) == grid.ny())) && (view_i(jnode) != 1)) {
          hmaskView(jnode, 0) = 0;
        } else {
          hmaskView(jnode, 0) = 1;
        }
      }
    }

    // Add field
    extraFields_->add(hmask);
  } else if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
    // Cubed-sphere grid
    atlas::functionspace::NodeColumns fs(functionSpace_);
    atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
      | atlas::option::levels(1));
    auto hmaskView = atlas::array::make_view<int, 2>(hmask);
    auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
    for (atlas::idx_t jnode = 0; jnode < hmask.shape(0); ++jnode) {
      hmaskView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
    }

    // Add field
    extraFields_->add(hmask);
  }

  // Print summary
  this->print(oops::Log::info());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_), halo_(other.halo_),
  grid_(other.grid_), unstructuredGrid_(other.unstructuredGrid_), partitioner_(other.partitioner_),
  mesh_(other.mesh_), levels_(other.levels_), lev2d_(other.lev2d_), vunit_(other.vunit_)  {
  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
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
  atlas::Field gmask = (*other.extraFields())["gmask"];
  extraFields_->add(gmask);
  gmaskSize_ = other.gmaskSize_;
}
// -------------------------------------------------------------------------------------------------
size_t Geometry::variableSize(const std::string & var) const {
  size_t levels = levels_;
  if (var.size() > 3) {
    if (var.substr(var.size()-3) == "_2d" || var.substr(var.size()-3) == "_2D") levels = 1;
  }
  return levels;
}
// -------------------------------------------------------------------------------------------------
size_t Geometry::maskLevel(const std::string & var, const size_t & level) const {
  size_t maskLevel = level;
  if (var.size() > 3) {
    if (var.substr(var.size()-3) == "_2d" || var.substr(var.size()-3) == "_2D") {
      if (lev2d_ == "first") {
        maskLevel = 1;
      } else if (lev2d_ == "last") {
        maskLevel = levels_-1;
      } else {
        ABORT("wrong lev2d value (first or last)");
      }
    }
  }
  return maskLevel;
}
// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes;
  for (const auto & var : vars.variables()) {
    sizes.push_back(this->variableSize(var));
  }
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix <<  "Quench geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (!unstructuredGrid_) {
    os << prefix << "Partitioner:" << std::endl;
    os << prefix << "- type: " << partitioner_.type() << std::endl;
  }
  os << prefix << "Function space:" << std::endl;
  os << prefix << "- type: " << functionSpace_.type() << std::endl;
  os << prefix << "- halo: " << halo_ << std::endl;
  os << prefix << "Vertical levels: " << std::endl;
  os << prefix << "- number: " << levels_ << std::endl;
  os << prefix << "- vunit: " << vunit_ << std::endl;
  os << prefix << "Mask size: " << static_cast<int>(gmaskSize_*100.0) << "%" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace quench
