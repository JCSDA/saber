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

  // Deal with poles for structured grids
  if (((gridType_ == "structured") || (gridType_ == "regular_lonlat")) && grid_.domain().global()) {
    // Get structured function space and grid
    const atlas::functionspace::StructuredColumns fs(functionSpace_);
    const atlas::StructuredGrid & grid = fs.grid();

    // Check whether the grid has duplicated points at the poles
    bool duplicated = false;
    int j = 0;
    if (grid.nx(j) > 1) {
      double lonlatPoint[] = {0, 0};
      grid.lonlat(0, j, lonlatPoint);
      double latRef = lonlatPoint[1];
      for (int i = 1; i < grid.nx(j); ++i) {
        grid.lonlat(i, j, lonlatPoint);
        duplicated = duplicated || (lonlatPoint[1] == latRef);
      }
    }
    j = grid.ny()-1;
    if (grid.nx(j) > 1) {
      double lonlatPoint[] = {0, 0};
      grid.lonlat(0, j, lonlatPoint);
      double latRef = lonlatPoint[1];
      for (int i = 1; i < grid.nx(j); ++i) {
        grid.lonlat(i, j, lonlatPoint);
        duplicated = duplicated || (lonlatPoint[1] == latRef);
      }
    }

    // Only keep the first of duplicated points in the owned mask
    const auto view_i = atlas::array::make_indexview<int, 1>(fs.index_i());
    const auto view_j = atlas::array::make_indexview<int, 1>(fs.index_j());
    auto ownedView = atlas::array::make_view<int, 2>(fieldsetOwnedMask["owned"]);
    for (int jnode = 0; jnode < fs.size(); ++jnode) {
      if (((view_j(jnode) == 0) || (view_j(jnode) == grid.ny()-1)) && (view_i(jnode) > 0)) {
        ownedView(jnode, 0) = 0;
      }
    }
  }

  // Setup geometry fields
  fields_ = atlas::FieldSet();

  // Add owned points mask -- this mask does not depend on the group so was precomputed
  fields_.add(fieldsetOwnedMask["owned"]);

  // Levels direction
  levelsAreTopDown_ = params.levelsAreTopDown.value();

  // Levels counter origin
  levelsCountFrom_ = params.levelsCountFrom.value();

  // Model data
  modelData_ = params.modelData.value();

  // Variable name alias
  setupAlias(params);

  // IO parameters
  io_ = params.io.value();

  // Interpolation
  const auto &interpParams = params.interpolation.value();
  if (interpParams != boost::none) {
    interpolation_ = interpParams->toConfiguration();
  } else {
    interpolation_ = eckit::LocalConfiguration();
    interpolation_.set("interpolation type", "unstructured");
  }

  // Check for duplicate points
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
  const auto ownedView = atlas::array::make_view<int, 2>(fields_["owned"]);
  size_t duplicatedPointsCount = 0;
  for (atlas::idx_t jnode = 0; jnode < fields_["owned"].shape(0); ++jnode) {
    // Duplicate point = owned==0 and ghost==0 (see util::setupFunctionSpace in oops)
    if (ghostView(jnode) == 0 && ownedView(jnode, 0) == 0) {
      ++duplicatedPointsCount;
    }
  }
  comm_.allReduceInPlace(duplicatedPointsCount, eckit::mpi::sum());
  duplicatePoints_ = (duplicatedPointsCount > 0);

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : params.groups.value()) {
    // Define group
    groupData group;

    // Copy group parameters
    group.params_ = groupParams;

    // Copy group index
    group.index_ = groupIndex;

    // Use this group index for all the group variables
    for (const auto & var : groupParams.variables.value()) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        throw eckit::UserError(
          "Same variable present in distinct groups " + var, Here());
      } else {
        groupIndex_[var] = groupIndex;
      }
    }

    // Number of levels
    group.levels_ = groupParams.levels.value();

    // Save group
    groups_.push_back(group);

    // Increment group index
    groupIndex++;
  }

  // Vertical coordinate
  for (auto & group : groups_) {
    setupVertCoord(group);
  }

  // Setup mask
  for (auto & group : groups_) {
    setupMask(group);
  }

  // Check lon/lat from files
  const auto &checkLonLatConf = params.checkLonLat.value();
  if (checkLonLatConf != boost::none) {
    checkLonLat(*checkLonLatConf);
  }

  // Print summary
  print(oops::Log::info());

  oops::Log::trace() << classname() << "::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t Geometry::groupIndex(const std::string & var) const {
  oops::Log::trace() << classname() << "::groupIndex starting" << std::endl;

  if (groupIndex_.find(var) == groupIndex_.end()) {
    throw eckit::Exception("cannot find group index for variable " + var, Here());
  }
  const size_t groupIndex = groupIndex_.at(var);

  oops::Log::trace() << classname() << "::groupIndex done" << std::endl;
  return groupIndex;
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

std::vector<size_t> Geometry::variableSizes(const std::vector<std::string> & varNames) const {
  oops::Log::trace() << classname() << "::variableSizes starting" << std::endl;

  // Create variables
  const oops::Variables vars(varNames);

  oops::Log::trace() << classname() << "::variableSizes done" << std::endl;
  return variableSizes(vars);
}

// -----------------------------------------------------------------------------

Interpolation & Geometry::getInterpolation(const Geometry & tgtGeom) const {
  oops::Log::trace() << classname() << "::getInterpolation starting" << std::endl;

  // Get target geometry UID (grid + "_" + paritioner + "@" + interpolation type)
  const std::string interpolationType = tgtGeom.interpolation().getString("interpolation type");
  const std::string tgtGeomUid = tgtGeom.grid().uid() + "_" + tgtGeom.partitioner().type()
    + "@" + interpolationType;

  // Look for this UID in the existing interpolations
  const auto it = interpolations_.find(tgtGeomUid);

  if (it != interpolations_.end()) {
    // Found existing interpolation
    oops::Log::info() << "Info     : Found existing interpolation for UID: " << tgtGeomUid
      << std::endl;

    // Return interpolation
    oops::Log::trace() << classname() << "::getInterpolation done" << std::endl;
    return *(it->second);
  } else {
    // Create GeometryData if needed
    if ((interpolationType == "unstructured") && !geomData_) {
      geomData_.reset(new oops::GeometryData(functionSpace_, fields_, levelsAreTopDown_, comm_));
    }

    // Create new interpolation
    std::shared_ptr<Interpolation> interpolation(new Interpolation(*this, tgtGeom));

    // Print interpolation type
    oops::Log::info() << "Info     : New interpolation created for UID: " << tgtGeomUid
      << std::endl;

    // Store interpolation
    interpolations_.insert({tgtGeomUid, interpolation});

    // Return interpolation
    oops::Log::trace() << classname() << "::getInterpolation done" << std::endl;
    return *interpolation;
  }
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
  if (duplicatePoints_) {
    os << prefix << "Duplicated points detected" << std::endl;
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
    os << prefix << "  - vertCoord: " << groups_[groupIndex].vertCoordAvg_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }

  if (!modelData_.empty()) {
    os << prefix << "Model data: " << modelData_ << std::endl;
  }

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::setupAlias(const GeometryParameters & params) {
  oops::Log::trace() << classname() << "::setupAlias starting" << std::endl;

  for (const auto & item : params.alias.value()) {
    eckit::LocalConfiguration confItem;
    item.serialize(confItem);
    alias_.push_back(confItem);
  }

  // Check alias consistency
  std::vector<std::string> vars;
  for (const auto & groupParams : params.groups.value()) {
    const std::vector<std::string> grpVars = groupParams.variables.value();
    vars.insert(vars.end(), grpVars.begin(), grpVars.end());
  }
  for (const auto & item : alias_) {
    const std::string codeVar = item.getString("in code");
    if (std::find(vars.begin(), vars.end(), codeVar) == vars.end()) {
      // Code variable not available in the list of variables anymore
      throw eckit::UserError("Alias error: code variable not available anymore", Here());
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

  oops::Log::trace() << classname() << "::setupAlias done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::setupVertCoord(groupData & group) {
  oops::Log::trace() << classname() << "::setupVertCoord starting" << std::endl;

  // Get optional parameters
  const auto &vertCoordConf = group.params_.vertCoordConf.value();

  // Get vertical coordinate name
  std::string vertCoordName = "vert_coord_" + std::to_string(group.index_);
  if (vertCoordConf != boost::none) {
    if (vertCoordConf->has("name")) {
      vertCoordName = vertCoordConf->getString("name");
    }
  }

  // Create vertical coordinate
  group.vertCoord_ = functionSpace_.createField<double>(
    atlas::option::name(vertCoordName) | atlas::option::levels(group.levels_));

  // Set interpolation metadata
  group.vertCoord_.metadata().set("interp_type", "default");

  // Get view
  auto vertCoordView = atlas::array::make_view<double, 2>(group.vertCoord_);

  if (vertCoordConf != boost::none) {
    // Vertical coordinate from a configuration
    if (vertCoordConf->has("profile")) {
      // From a vector of doubles (one for each level)
      const std::vector<double> profile = vertCoordConf->getDoubleVector("profile");
      if (profile.size() != group.levels_) {
        throw eckit::UserError("Wrong number of levels in the user-specified vertical coordinate",
          Here());
      }
      for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vertCoordView(jnode, jlevel) = profile[jlevel];
        }
      }
    } else {
      // Get variable to read
      const std::string varName = vertCoordConf->getString("variable");
      const oops::Variables vertCoordVars(std::vector<std::string>({varName}));

      // Add group index for this variable
      groupIndex_[varName] = group.index_;

      // Create field
      Fields field(*this, vertCoordVars, util::DateTime(), false);

      // Read field
      field.read(*vertCoordConf);

      // Get view
      const auto view = atlas::array::make_view<double, 2>(field.fieldSet()[varName]);

      // Copy 3D field
      for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
        for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
          vertCoordView(jnode, jlevel) = view(jnode, jlevel);
        }
      }
    }
  } else {
    // From level index (default)
    for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        vertCoordView(jnode, jlevel) = static_cast<double>(jlevel+levelsCountFrom_);
      }
    }
  }

  // Get owned view
  const auto ownedView = atlas::array::make_view<int, 2>(fields_["owned"]);

  // Average vertical coordinate
  for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
    // Initialization
    double avg = 0.0;
    double counter = 0.0;

    // Loop over owned points
    for (atlas::idx_t jnode = 0; jnode < group.vertCoord_.shape(0); ++jnode) {
      if (ownedView(jnode, 0) == 1) {
        avg += vertCoordView(jnode, jlevel);
        counter += 1.0;
      }
    }

    // Communication
    comm_.allReduceInPlace(avg, eckit::mpi::sum());
    comm_.allReduceInPlace(counter, eckit::mpi::sum());

    // Normalization
    if (counter > 0.0) {
      avg /= counter;
    }

    // Update profile
    group.vertCoordAvg_.push_back(avg);
  }

  // Add orography (mountain) on bottom level
  const auto &orographyParams = group.params_.orography.value();
  if (orographyParams != boost::none) {
    // Get top latitude value
    const atlas::PointLonLat topPoint({orographyParams->topLon.value(),
      orographyParams->topLat.value()});

    // Get lon/lat view
    const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());

    for (atlas::idx_t jnode = 0; jnode < lonlatView.shape(0); ++jnode) {
      // Get delta
      const double delta = (group.levels_ == 1) ? 1.0 :
        vertCoordView(jnode, group.levels_-2)-vertCoordView(jnode, group.levels_-1);

      // Get x and y points
      const atlas::PointLonLat xPoint({lonlatView(jnode, 0), orographyParams->topLat.value()});
      const atlas::PointLonLat yPoint({orographyParams->topLon.value(), lonlatView(jnode, 1)});

      // Compute normalization
      const double dxNorm = atlas::util::Earth().distance(xPoint, topPoint)
        /orographyParams->zonalLength.value();
      const double dyNorm = atlas::util::Earth().distance(yPoint, topPoint)
        /orographyParams->meridionalLength.value();
      const double distNorm = std::sqrt(dxNorm*dxNorm+dyNorm*dyNorm);

      // Define orography
      const double orography = delta*orographyParams->height.value()*oops::gc99(distNorm);

      // Add orography to existing vertical coordinate
      vertCoordView(jnode, group.levels_-1) += orography;
    }
  }

  // Add vertical coordinate in Geometry fields
  fields_.add(group.vertCoord_);

  oops::Log::trace() << classname() << "::setupVertCoord starting" << std::endl;
}


// -----------------------------------------------------------------------------

void Geometry::setupMask(groupData & group) {
  oops::Log::trace() << classname() << "::setupMask starting" << std::endl;

  // Default mask, set to 1 (true)
  const std::string gmaskName = "gmask_" + std::to_string(group.index_);
  atlas::Field gmask = functionSpace_.createField<int>(
    atlas::option::name(gmaskName) | atlas::option::levels(group.levels_));
  auto maskView = atlas::array::make_view<int, 2>(gmask);
  maskView.assign(1);

  // Owned view
  auto ownedView = atlas::array::make_view<int, 2>(fields_["owned"]);

  // Specific mask
  if (group.params_.maskType.value() == "none") {
    // No mask
  } else if (group.params_.maskType.value() == "sea") {
    // Read sea mask

    // Lon/lat sizes
    size_t nlon = 0;
    size_t nlat = 0;

    // File path
    ASSERT(group.params_.maskPath.value());
    const std::string ncFilePath = *group.params_.maskPath.value();

    // NetCDF IDs
    int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

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
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Broadcast coordinates and land-sea mask
    comm_.broadcast(lon.begin(), lon.end(), 0);
    comm_.broadcast(lat.begin(), lat.end(), 0);
    comm_.broadcast(lsm.begin(), lsm.end(), 0);

    // Build KD-tree
    const atlas::Geometry geometry(atlas::util::Earth::radius());
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
      const atlas::functionspace::StructuredColumns fs(functionSpace_);
      const auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
      auto maskView = atlas::array::make_view<int, 2>(gmask);
      for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
        if (ownedView(jnode, 0) == 1) {
          // Find nearest neighbor
          size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
            lonlatView(jnode, 1)}).payload();

          // Ocean points for all levels
          for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
            if (lsm[nn] == 0) {
              maskView(jnode, jlevel) = 1;
            } else {
              maskView(jnode, jlevel) = 0;
            }
          }
        }
      }
    } else {
      throw eckit::NotImplemented("Sea mask not supported for " + functionSpace_.type() + " yet",
        Here());
    }
  } else {
    throw eckit::UserError("Wrong mask type", Here());
  }

  // Mask size
  group.gmaskSize_ = 0.0;
  size_t domainSize = 0.0;
  for (atlas::idx_t jnode = 0; jnode < gmask.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < gmask.shape(1); ++jlevel) {
      if (ownedView(jnode, 0) == 1) {
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

  // Add mask in Geometry fields
  fields_.add(gmask);

  oops::Log::trace() << classname() << "::setupMask done" << std::endl;
}

// -----------------------------------------------------------------------------

void Geometry::checkLonLat(const eckit::Configuration & checkLonLatConf) {
  oops::Log::trace() << classname() << "::checkLonLat starting" << std::endl;

  // Return if configuration is empty
  if (checkLonLatConf.empty()) {
    return;
  }

  // Get variable to read
  const std::string lonName = checkLonLatConf.getString("longitude", "longitude");
  const std::string latName = checkLonLatConf.getString("latitude", "latitude");
  const oops::Variables lonLatVars(std::vector<std::string>({lonName, latName}));

  // Add new group to read coordinates
  groupData coordGroup;
  coordGroup.levels_ = 1;
  groupIndex_["longitude"] = groups_.size();
  groupIndex_["latitude"] = groups_.size();
  groups_.push_back(coordGroup);

  // Create field
  Fields field(*this, lonLatVars, util::DateTime(), false);

  // Read field
  field.read(checkLonLatConf);

  // Get views
  const auto lonView = atlas::array::make_view<double, 2>(field.fieldSet()[lonName]);
  const auto latView = atlas::array::make_view<double, 2>(field.fieldSet()[latName]);

  // Get lon/lat view
  const auto lonlatView = atlas::array::make_view<double, 2>(functionSpace_.lonlat());

  // Get ghost view
  const auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // Get lon/lat check tolerance
  const double lonLatTol = checkLonLatConf.getDouble("tolerance", 1.0e-12);

  // Check lon/lat
  for (atlas::idx_t jnode = 0; jnode < functionSpace_.lonlat().shape(0); ++jnode) {
    if (ghostView(jnode) == 0) {
      if (std::abs(lonView(jnode, 0)-lonlatView(jnode, 0)) > lonLatTol) {
        throw eckit::Exception("inaccurate longitude", Here());
      }
      if (std::abs(latView(jnode, 0)-lonlatView(jnode, 1)) > lonLatTol) {
        throw eckit::Exception("inaccurate latitude", Here());
      }
    }
  }

  oops::Log::trace() << classname() << "::checkLonLat starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
