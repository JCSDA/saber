/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
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

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

#include "src/Fields.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); throw eckit::Exception(s + ": " + msg, Here());}

namespace quench {

// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), groups_() {
  oops::Log::trace() << classname() << "::Geometry starting" << std::endl;

  GeometryParameters params;
  params.deserialize(config);

  // Setup atlas geometric data structures
  atlas::FieldSet fieldsetOwnedMask;
  util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh_, functionSpace_,
    fieldsetOwnedMask);
  halo_ = params.halo.value();
  gridType_ = params.grid.value().getString("type", "no_type");

  // Setup geometry fields
  fields_ = atlas::FieldSet();

  // Add owned points mask -- this mask does not depend on the group so was precomputed
  fields_->add(fieldsetOwnedMask.field("owned"));

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : params.groups.value()) {
    // Use this group index for all the group variables
    for (const auto & var : groupParams.variables.value()) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        throw eckit::UserError(
          "Same variable present in distinct groups " + var, Here());
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

    // Vertical coordinate
    const boost::optional<std::vector<double>> &vert_coordParams = groupParams.vert_coord.value();
    const boost::optional<eckit::LocalConfiguration> &vert_coordParamsFromFile =
      groupParams.vert_coordFromFile.value();
    const std::string vert_coordName = "vert_coord_" + std::to_string(groupIndex);
    group.vert_coord_ = functionSpace_.createField<double>(
      atlas::option::name(vert_coordName) | atlas::option::levels(group.levels_));
    group.vert_coord_.metadata().set("interp_type", "default");
    auto vert_coordView = atlas::array::make_view<double, 2>(group.vert_coord_);
    if (vert_coordParams != boost::none) {
      // From a vector of doubles (one for each level)
      if (vert_coordParams->size() != group.levels_) {
        throw eckit::UserError("Wrong number of levels in the user-specified vertical coordinate",
          Here());
      }
      for (atlas::idx_t jnode = 0; jnode < group.vert_coord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vert_coordView(jnode, jlevel) = (*vert_coordParams)[jlevel];
        }
      }
    } else if (vert_coordParamsFromFile != boost::none) {
      // From a file
      const std::vector<std::string> vert_coordVars =
        vert_coordParamsFromFile->getStringVector("variables");
      const oops::Variables vert_coordVar(vert_coordVars);
      eckit::LocalConfiguration fileGeomConfig(config);
      std::vector<eckit::LocalConfiguration> groupsConfig(1);
      groupsConfig[0].set("variables", vert_coordVars);
      groupsConfig[0].set("levels", 1);
      fileGeomConfig.set("groups", groupsConfig);
      Geometry fileGeom(fileGeomConfig);
      Fields field(fileGeom, vert_coordVar, util::DateTime());
      field.read(*vert_coordParamsFromFile);
      const auto view = atlas::array::make_view<double, 2>(field.fieldSet()[vert_coordVars[0]]);
      for (atlas::idx_t jnode = 0; jnode < group.vert_coord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vert_coordView(jnode, jlevel) = view(jnode, jlevel);
        }
      }
    } else {
      // From level index
      for (atlas::idx_t jnode = 0; jnode < group.vert_coord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vert_coordView(jnode, jlevel) = static_cast<double>(jlevel+1);
        }
      }
    }

    // Average vertical coordinate
    const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
    const auto ownedView = atlas::array::make_view<int, 2>(fields_.field("owned"));
    for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
      double avg = 0.0;
      double counter = 0.0;
      for (atlas::idx_t jnode = 0; jnode < group.vert_coord_.shape(0); ++jnode) {
        if (ghostView(jnode) == 0 && ownedView(jnode, 0) == 1) {
          avg += vert_coordView(jnode, jlevel);
          counter += 1.0;
        }
      }
      comm_.allReduceInPlace(avg, eckit::mpi::sum());
      comm_.allReduceInPlace(counter, eckit::mpi::sum());
      if (counter > 0.0) {
        avg /= counter;
      }
      group.vert_coord_avg_.push_back(avg);
    }

    // Add orography (mountain) on bottom level
    const boost::optional<OrographyParameters> &orographyParams = groupParams.orography.value();
    if (orographyParams != boost::none) {
      const atlas::PointLonLat topPoint({orographyParams->topLon.value(),
        orographyParams->topLat.value()});
      const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
      for (atlas::idx_t jnode = 0; jnode < lonlatView.shape(0); ++jnode) {
        const double delta = (group.levels_ == 1) ? 1.0 :
          vert_coordView(jnode, group.levels_-2)-vert_coordView(jnode, group.levels_-1);
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
    fields_->add(group.vert_coord_);

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
      ASSERT(groupParams.maskPath.value());
      readSeaMask(*groupParams.maskPath.value(), group.levels_, group.lev2d_, gmask);
    } else {
      throw eckit::UserError("Wrong mask type", Here());
    }
    fields_->add(gmask);

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

  // Levels direction
  levelsAreTopDown_ = params.levelsAreTopDown.value();

  // Model data
  modelData_ = params.modelData.value();

  // Alias
  for (const auto & item : params.alias.value()) {
    eckit::LocalConfiguration confItem;
    item.serialize(confItem);
    alias_.push_back(confItem);
  }

  // Check alias consistency
  std::vector<std::string> vars;
  for (const auto & groupParams : params.groups.value()) {
    const std::vector grpVars = groupParams.variables.value();
    vars.insert(vars.end(), grpVars.begin(), grpVars.end());
  }
  for (const auto & item : alias_) {
    const std::string codeVar = item.getString("in code");
    if (std::find(vars.begin(), vars.end(), codeVar) == vars.end()) {
      // Code variable not available in the list of variables anymore
      throw eckit::UserError("Alias error: duplicated code variable", Here());
    } else {
      // Remove code variable from the list of available variables
      vars.erase(std::remove(vars.begin(), vars.end(), codeVar), vars.end());
    }
  }
  for (const auto & item : alias_) {
    const std::string fileVar = item.getString("in file");
    if (std::find(vars.begin(), vars.end(), fileVar) == vars.end()) {
      // Add file variable to the list of variables
      vars.push_back(fileVar);
    } else {
      // File variable is already present in the list of variables
      throw eckit::UserError("Alias error: duplicated file variable", Here());
    }
  }

  // Latitudes from south to north in files
  latSouthToNorth_ = params.latSouthToNorth.value();

  // Interpolation
  const boost::optional<InterpolationParameters> &interpParams = params.interpolation.value();
  if (interpParams != boost::none) {
    interpolation_ = interpParams->toConfiguration();
  } else {
    interpolation_ = eckit::LocalConfiguration();
      interpolation_.set("interpolation type", "unstructured");
  }

  // Check for duplicate points
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
  const auto ownedView = atlas::array::make_view<int, 2>(fields_.field("owned"));
  size_t duplicatedPointsCount = 0;
  for (atlas::idx_t jnode = 0; jnode < fields_.field("owned").shape(0); ++jnode) {
    // Duplicate point = owned==0 and ghost==0 (see util::setupFunctionSpace in oops)
    if (ghostView(jnode) == 0 && ownedView(jnode, 0) == 0) {
      ++duplicatedPointsCount;
    }
  }
  comm_.allReduceInPlace(duplicatedPointsCount, eckit::mpi::sum());
  duplicatePoints_ = (duplicatedPointsCount > 0);

  // GeometryData
  if (interpolation_.getString("interpolation type") == "unstructured") {
    geomData_.reset(new oops::GeometryData(functionSpace_, fields_, levelsAreTopDown_, comm_));
  }

  // Print summary
  this->print(oops::Log::info());

  oops::Log::trace() << classname() << "::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other)
  : comm_(other.comm_), halo_(other.halo_), grid_(other.grid_), gridType_(other.gridType_),
  partitioner_(other.partitioner_), mesh_(other.mesh_), groupIndex_(other.groupIndex_),
  levelsAreTopDown_(other.levelsAreTopDown_), modelData_(other.modelData_), alias_(other.alias_),
  latSouthToNorth_(other.latSouthToNorth_), interpolation_(other.interpolation_),
  duplicatePoints_(other.duplicatePoints_) {
  oops::Log::trace() << classname() << "::Geometry starting" << std::endl;

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
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported",
      Here());
  } else {
    throw eckit::NotImplemented(other.functionSpace_.type() + " function space not supported yet",
      Here());
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

    // Copy averaged vertical coordinate
    group.vert_coord_avg_ = other.groups_[groupIndex].vert_coord_avg_;

    // Copy mask size
    group.gmaskSize_ = other.groups_[groupIndex].gmaskSize_;

    // Save group
    groups_.push_back(group);
  }

  // Geometry data
  if (interpolation_.getString("interpolation type") == "unstructured") {
    geomData_.reset(new oops::GeometryData(functionSpace_, fields_, levelsAreTopDown_, comm_));
  }

  oops::Log::trace() << classname() << "::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::variableSizes starting" << std::endl;

  std::vector<size_t> sizes;
  for (const auto & var : vars) {
    sizes.push_back(levels(var.name()));
  }

  oops::Log::trace() << classname() << "::variableSizes done" << std::endl;
  return sizes;
}

// -----------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

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
    os << prefix << "  - vert_coord: " << groups_[groupIndex].vert_coord_avg_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }

  if (!modelData_.empty()) {
    os << prefix << "Model data: " << modelData_ << std::endl;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::readSeaMask(const std::string & maskPath,
                           const size_t & levels,
                           const std::string & lev2d,
                           atlas::Field & gmask) const {
  oops::Log::trace() << classname() << "::readSeaMask starting" << std::endl;

  // Lon/lat sizes
  size_t nlon = 0;
  size_t nlat = 0;

  // NetCDF IDs
  int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(maskPath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, maskPath);

    // Get lon/lat sizes
    if ((retval = nc_inq_dimid(ncid, "lon", &nlon_id))) ERR(retval, "lon");
    if ((retval = nc_inq_dimid(ncid, "lat", &nlat_id))) ERR(retval, "lat");
    if ((retval = nc_inq_dimlen(ncid, nlon_id, &nlon))) ERR(retval, "lon");
    if ((retval = nc_inq_dimlen(ncid, nlat_id, &nlat))) ERR(retval, "lat");
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
    if ((retval = nc_inq_varid(ncid, "lon", &lon_id))) ERR(retval, "lon");
    if ((retval = nc_inq_varid(ncid, "lat", &lat_id))) ERR(retval, "lat");
    if ((retval = nc_inq_varid(ncid, "LSMASK", &lsm_id))) ERR(retval, "LMASK");

    // Read data
    std::vector<float> zlon(nlon);
    std::vector<float> zlat(nlat);
    std::vector<uint8_t> zlsm(nlat*nlon);
    if ((retval = nc_get_var_float(ncid, lon_id, zlon.data()))) ERR(retval, "lon");
    if ((retval = nc_get_var_float(ncid, lat_id, zlat.data()))) ERR(retval, "lat");
    if ((retval = nc_get_var_ubyte(ncid, lsm_id, zlsm.data()))) ERR(retval, "LMASK");

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
    if ((retval = nc_close(ncid))) ERR(retval, maskPath);
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
    throw eckit::NotImplemented("Sea mask not supported for " + functionSpace_.type() + " yet",
      Here());
  }

  oops::Log::trace() << classname() << "::readSeaMask done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
