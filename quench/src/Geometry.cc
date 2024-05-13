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

#include "oops/generic/gc99.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & config, const eckit::mpi::Comm & comm)
  : comm_(comm), groups_()
{
  GeometryParameters params;
  params.deserialize(config);

  // Setup atlas geometric data structures
  atlas::FieldSet fieldsetOwnedMask;
  util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh_,
                           functionSpace_, fieldsetOwnedMask);

  halo_ = params.halo.value();
  gridType_ = params.grid.value().getString("type", "no_type");

  // Setup geometry fields
  fields_ = atlas::FieldSet();

  // Add owned points mask -- this mask does not depend on the group so was precomputed
  fields_->add(fieldsetOwnedMask.field("owned"));

  if (!grid_.domain().global()) {
    // Area
    atlas::Field area = functionSpace_.createField<double>(
      atlas::option::name("area") | atlas::option::levels(1));
    auto areaView = atlas::array::make_view<double, 2>(area);
    const atlas::StructuredGrid grid(grid_);
    const atlas::functionspace::StructuredColumns fs(functionSpace_);
    const auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
    const auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
    for (atlas::idx_t jnode = 0; jnode < area.shape(0); ++jnode) {
      // Initialization
      int i = view_i(jnode)-1;
      int j = view_j(jnode)-1;
      double dist_i = 0.0;
      double dist_j = 0.0;

      // i-direction component
      if (i == 0) {
        dist_i = atlas::util::Earth().distance(grid.lonlat(i, j), grid.lonlat(i+1, j));
      } else if (i == grid.nx(j)-1) {
        dist_i = atlas::util::Earth().distance(grid.lonlat(i-1, j), grid.lonlat(i, j));
      } else {
        dist_i = 0.5*atlas::util::Earth().distance(grid.lonlat(i-1, j), grid.lonlat(i+1, j));
      }

      // j-direction component
      if (j == 0) {
        dist_j = atlas::util::Earth().distance(grid.lonlat(i, j), grid.lonlat(i, j+1));
      } else if (j == grid.ny()-1) {
        dist_j = atlas::util::Earth().distance(grid.lonlat(i, j-1), grid.lonlat(i, j));
      } else {
        dist_j = 0.5*atlas::util::Earth().distance(grid.lonlat(i, j-1), grid.lonlat(i, j+1));
     }

      // Local scale
      areaView(jnode, 0) = dist_i*dist_j;
    }
    fields_->add(area);
  }

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

    // Average vertical coordinate
    const boost::optional<std::vector<double>> &vert_coordParams = groupParams.vert_coord.value();
    if (vert_coordParams != boost::none) {
      if (vert_coordParams->size() != group.levels_) {
        ABORT("Wrong number of levels in the user-specified vertical coordinate");
      }
    }
    for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
      if (vert_coordParams != boost::none) {
        group.vert_coord_.push_back((*vert_coordParams)[jlevel]);
      } else {
        group.vert_coord_.push_back(static_cast<double>(jlevel+1));
      }
    }

    // Vertical coordinate field
    const std::string vert_coordName = "vert_coord_" + std::to_string(groupIndex);
    atlas::Field vert_coord = functionSpace_.createField<double>(
      atlas::option::name(vert_coordName) | atlas::option::levels(group.levels_));
    auto vert_coordView = atlas::array::make_view<double, 2>(vert_coord);
    for (atlas::idx_t jnode = 0; jnode < vert_coord.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        vert_coordView(jnode, jlevel) = group.vert_coord_[jlevel];
      }
    }

    // Add orography (mountain) on bottom level
    const boost::optional<OrographyParameters> &orographyParams = groupParams.orography.value();
    if (orographyParams != boost::none) {
      const atlas::PointLonLat topPoint({orographyParams->topLon.value(),
        orographyParams->topLat.value()});
      const double delta = (group.levels_ == 1) ? 1.0 :
        group.vert_coord_[group.levels_-2]-group.vert_coord_[group.levels_-1];
      const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
      for (atlas::idx_t jnode = 0; jnode < lonlatView.shape(0); ++jnode) {
        const atlas::PointLonLat xPoint({lonlatView(jnode, 0), orographyParams->topLat.value()});
        const atlas::PointLonLat yPoint({orographyParams->topLon.value(), lonlatView(jnode, 1)});
        double dxNorm = atlas::util::Earth().distance(xPoint, topPoint)
          /orographyParams->zonalLength.value();
        double dyNorm = atlas::util::Earth().distance(yPoint, topPoint)
          /orographyParams->meridionalLength.value();
        double distNorm = std::sqrt(dxNorm*dxNorm+dyNorm*dyNorm);
        double orography = delta*orographyParams->height.value()*oops::gc99(distNorm);
        vert_coordView(jnode, group.levels_-1) += orography;
      }
    }
    fields_->add(vert_coord);

    // Default mask, set to 1 (true)
    const std::string gmaskName = "gmask_" + std::to_string(groupIndex);
    atlas::Field gmask = functionSpace_.createField<int>(
      atlas::option::name(gmaskName) | atlas::option::levels(group.levels_));
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
    fields_->add(gmask);

    // Mask size
    group.gmaskSize_ = 0.0;
    size_t domainSize = 0.0;
    auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
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
  grid_(other.grid_), gridType_(other.gridType_), partitioner_(other.partitioner_),
  mesh_(other.mesh_), groupIndex_(other.groupIndex_)  {
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
    ABORT(other.functionSpace_.type() + " function space not supported");
  } else {
    ABORT(other.functionSpace_.type() + " function space not supported yet");
  }

  // Copy geometry fields
  fields_ = util::shareFields(other.fields_);

  // Copy groups
  for (size_t groupIndex = 0; groupIndex < other.groups_.size(); ++groupIndex) {
    // Define group
    groupData group;

    // Copy number of levels
    group.levels_ = other.groups_[groupIndex].levels_;

    // Copy corresponding level for 2D variables (first or last)
    group.lev2d_ = other.groups_[groupIndex].lev2d_;

    // Copy vertical coordinate
    group.vert_coord_ = other.groups_[groupIndex].vert_coord_;

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
void Geometry::print(std::ostream & os) const {
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix <<  "Quench geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (!grid_.domain().global()) {
    os << prefix << "Regional grid detected" << std::endl;
  }
  if (partitioner_) {
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
    os << prefix << "  - vert_coord: " << groups_[groupIndex].vert_coord_ << std::endl;
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
    std::vector<float> zlon(nlon);
    std::vector<float> zlat(nlat);
    std::vector<uint8_t> zlsm(nlat*nlon);
    if ((retval = nc_get_var_float(ncid, lon_id, zlon.data()))) ERR(retval);
    if ((retval = nc_get_var_float(ncid, lat_id, zlat.data()))) ERR(retval);
    if ((retval = nc_get_var_ubyte(ncid, lsm_id, zlsm.data()))) ERR(retval);

    // Copy data
    for (size_t ilon = 0; ilon < nlon; ++ilon) {
      lon[ilon] = static_cast<double>(zlon[ilon]);
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
      lat[ilat] = static_cast<double>(zlat[ilat]);
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
     for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lsm[ilat*nlon+ilon] = static_cast<int>(zlsm[ilat*nlon+ilon]);
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
