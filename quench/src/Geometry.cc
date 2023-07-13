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

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Geometry::Geometry(const Parameters_ & params,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), gridType_("no_type"), groups_() {
  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Halo
  const boost::optional<size_t> &halo = params.halo.value();
  if (halo != boost::none) {
    halo_ = *halo;
  } else {
    halo_ = 0;
  }

  // Set flag
  unstructuredGrid_ = false;

  // Setup grid
  eckit::LocalConfiguration gridParams = params.grid.value();
  oops::Log::info() << "Info     : Grid config: " << gridParams << std::endl;
  if (gridParams.has("type")) {
    gridType_ = gridParams.getString("type");
    if (gridType_ == "unstructured") {
      // Split unstructured grid among processors
      std::vector<double> xyFull = gridParams.getDoubleVector("xy");
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
        if (params.noPointOnLastTask.value() && (comm_.size() > 1)) {
          if (rank == comm_.size()-1) rank = 0;
        } else {
          if (rank == comm_.size()) rank = 0;
        }
      }

      // Reset coordinates
      gridParams.set("xy", xy);

      // Set flag
      unstructuredGrid_ = true;
    }
  }
  grid_ = atlas::Grid(gridParams);

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

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : params.groups.value()) {
    // Use this group index for all the group variables
    for (const auto & var : groupParams.variables.value()) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        ABORT("Same variable present in distinct groups");
      } else {
        groupIndex_[var] = groupIndex;
      }
    }

    // Define group
    groupData group;

    // Number of levels
    group.levels_ = groupParams.levels.value();

    // Corresponding level for 2D variables (first or last)
    group.lev2d_ = groupParams.lev2d.value();

    // Vertical unit
    const boost::optional<std::vector<double>> &vunitParams = groupParams.vunit.value();
    if (vunitParams != boost::none) {
      if (vunitParams->size() != group.levels_) {
        ABORT("Wrong number of levels in the user-specified vertical unit");
      }
    }
    for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
      if (vunitParams != boost::none) {
        group.vunit_.push_back((*vunitParams)[jlevel]);
      } else {
        group.vunit_.push_back(static_cast<double>(jlevel+1));
      }
    }

    // Default mask, set to 1 (true)
    atlas::Field gmask = functionSpace_.createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(group.levels_));
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    maskView.assign(1);

    // Specific mask
    if (groupParams.maskType.value() == "none") {
      // No mask
    } else if (groupParams.maskType.value() == "sea") {
      // Read sea mask
      readSeaMask(groupParams.maskPath.value(), group.levels_, group.lev2d_, gmask);
    } else {
      ABORT("Wrong mask type");
    }

    // Fill extra geometry fields
    group.extraFields_ = atlas::FieldSet();

    // Vertical unit
    atlas::Field vunit = functionSpace_.createField<double>(
      atlas::option::name("vunit") | atlas::option::levels(group.levels_));
    auto vunitView = atlas::array::make_view<double, 2>(vunit);
    for (atlas::idx_t jnode = 0; jnode < vunit.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
         vunitView(jnode, jlevel) = group.vunit_[jlevel];
      }
    }
    group.extraFields_->add(vunit);

    // Geographical mask
    group.extraFields_->add(gmask);

    // Halo mask
    if (functionSpace_.type() == "StructuredColumns") {
      // Structured columns
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      atlas::StructuredGrid grid = fs.grid();
      atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
        | atlas::option::levels(1));
      auto hmaskView = atlas::array::make_view<int, 2>(hmask);
      auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
      auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
          atlas::idx_t jnode = fs.index(i, j);
          hmaskView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }
      }

      // Special case for lon/lat grids
      if (grid_.name().compare(0, 1, std::string{"L"}) == 0) {
        for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
          for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (((view_j(jnode) == 1) || (view_j(jnode) == grid.ny())) && (view_i(jnode) != 1)) {
              hmaskView(jnode, 0) = 0;
            }
          }
        }
      }

      // Add halo mask
      group.extraFields_->add(hmask);
    } else if (functionSpace_.type() == "NodeColumns") {
      // NodeColumns
      if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::NodeColumns fs(functionSpace_);
        atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
          | atlas::option::levels(1));
        auto hmaskView = atlas::array::make_view<int, 2>(hmask);
        auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
        for (atlas::idx_t jnode = 0; jnode < hmask.shape(0); ++jnode) {
          hmaskView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }

        // Add halo mask
        group.extraFields_->add(hmask);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(functionSpace_);
        atlas::Field hmask = fs.createField<int>(atlas::option::name("hmask")
          | atlas::option::levels(1));
        auto hmaskView = atlas::array::make_view<int, 2>(hmask);
        auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
        for (atlas::idx_t jnode = 0; jnode < hmask.shape(0); ++jnode) {
          hmaskView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }

        // Add halo mask
        group.extraFields_->add(hmask);
      }
    }

    // Mask size
    group.gmaskSize_ = 0.0;
    size_t domainSize = 0.0;
    for (atlas::idx_t jnode = 0; jnode < gmask.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < gmask.shape(1); ++jlevel) {
        if (ghostView(jnode) == 0) {
          if (maskView(jnode, jlevel) == 1) {
            group.gmaskSize_ += 1.0;
          }
          domainSize++;
        }
      }
    }
    comm_.allReduceInPlace(group.gmaskSize_, eckit::mpi::sum());
    comm_.allReduceInPlace(domainSize, eckit::mpi::sum());
    if (domainSize > 0) {
      group.gmaskSize_ = group.gmaskSize_/static_cast<double>(domainSize);
    }

    // Save group
    groups_.push_back(group);

    // Increment group index
    groupIndex++;
  }

  // Print summary
  this->print(oops::Log::info());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_), halo_(other.halo_),
  grid_(other.grid_), gridType_(other.gridType_), unstructuredGrid_(other.unstructuredGrid_),
  partitioner_(other.partitioner_), mesh_(other.mesh_), groupIndex_(other.groupIndex_)  {
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

  // Copy groups
  for (size_t groupIndex = 0; groupIndex < other.groups_.size(); ++groupIndex) {
    // Define group
    groupData group;

    // Copy number of levels
    group.levels_ = other.groups_[groupIndex].levels_;

    // Copy corresponding level for 2D variables (first or last)
    group.lev2d_ = other.groups_[groupIndex].lev2d_;

    // Copy vertical unit
    group.vunit_ = other.groups_[groupIndex].vunit_;

    // Copy extra fields
    group.extraFields_ = atlas::FieldSet();
    group.extraFields_->add(other.groups_[groupIndex].extraFields_["vunit"]);
    group.extraFields_->add(other.groups_[groupIndex].extraFields_["gmask"]);
    if (other.groups_[groupIndex].extraFields_.has("hmask")) {
      group.extraFields_->add(other.groups_[groupIndex].extraFields_["hmask"]);
    }

    // Copy mask size
    group.gmaskSize_ = other.groups_[groupIndex].gmaskSize_;

    // Save group
    groups_.push_back(group);
  }
}
// -----------------------------------------------------------------------------
size_t Geometry::levels(const std::string & var) const {
  if (groupIndex_.count(var) == 0) {
    ABORT("Variable " + var + " not found in groupIndex_");
  }
  return groups_[groupIndex_.at(var)].levels_;
}
// -----------------------------------------------------------------------------
size_t Geometry::groupIndex(const std::string & var) const {
  if (groupIndex_.count(var) == 0) {
    ABORT("Variable " + var + " not found in groupIndex_");
  }
  return groupIndex_.at(var);
}
// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes;
  for (const auto & var : vars.variables()) {
    sizes.push_back(levels(var));
  }
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool includeHaloForRealLife) const {
  const auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  const auto ghost = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // TODO(Algo): Remove/fix the hack below when GeometryData local KD tree needs
  // to be set up correctly (e.g. when UnstructuredInterpolator is used).
  // For now never include halo in the latlon output because halo points from
  // some atlas grids (e.g. gaussian) can have unrealistic latitudes (e.g. more
  // than 90 degrees) and those latitudes can't be handled by KD trees.
  // Global KD trees created in GeometryData are used for communication and
  // don't need halo information.
  // Local KD trees in GeometryData need halo information but aren't used unless
  // UnstructuredInterpolator is used.
  bool includeHalo = false;
  const size_t npts = functionSpace_.size();
  const size_t nptsReturned = [&]() {
    if (includeHalo && comm_.size() > 1) {
      return npts;
    } else {
      size_t result = 0;
      for (atlas::idx_t i = 0; i < ghost.shape(0); ++i) {
        if (ghost(i) == 0) {
          result++;
        }
      }
      return result;
    }
  }();

  lats.resize(nptsReturned);
  lons.resize(nptsReturned);

  size_t count = 0;
  for (size_t jj = 0; jj < npts; ++jj) {
    // copy owned points, i.e. points with ghost==?
    if (ghost(jj) == 0 || (includeHalo && comm_.size() > 1)) {
      lats[count] = lonlat(jj, 1);
      lons[count] = lonlat(jj, 0);
      if (lons[count] < 0.0) lons[count] += 360.0;
      count++;
    }
  }
  ASSERT(count == nptsReturned);
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
  os << prefix << "Groups: " << std::endl;
  for (size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
    os << prefix << "- Group " << groupIndex << ":" << std::endl;
    os << prefix << "  Vertical levels: " << std::endl;
    os << prefix << "  - number: " << levels(groupIndex) << std::endl;
    os << prefix << "  - vunit: " << groups_[groupIndex].vunit_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }
}
// -----------------------------------------------------------------------------
void Geometry::readSeaMask(const std::string & maskPath,
                           const size_t & levels,
                           const std::string & lev2d,
                           atlas::Field & gmask) const {
  // Lon/lat sizes
  size_t nlon = 0;
  size_t nlat = 0;

  // NetCDF IDs
  int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(maskPath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

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

  // Ghost points
  atlas::Field ghost = functionSpace_.ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);

  if (functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(functionSpace_);
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Find nearest neighbor
        size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
          lonlatView(jnode, 1)}).payload();

        // Ocean points for all levels
        for (size_t jlevel = 0; jlevel < levels; ++jlevel) {
          if (lsm[nn] == 0) {
             maskView(jnode, jlevel) = 1;
           } else {
             maskView(jnode, jlevel) = 0;
           }
         }

        // Ocean + small islands for:
        // - the first level of 3D fields,
        // - the 2D fields if lev2d = "first"
        if (lsm[nn] == 3) {
          if ((levels > 1) || (lev2d == "first")) {
            maskView(jnode, 0) = 1;
          }
        }
      }
    }
  } else {
    ABORT("Sea mask not supported for " + functionSpace_.type() + " yet");
  }
}
// -----------------------------------------------------------------------------
}  // namespace quench
