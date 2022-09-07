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

#define ERR(e) {ABORT(nc_strerror(e));}

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
  const boost::optional<std::string> &iodaFile = params.iodaFile.value();
  if (gridParams != boost::none) {
    oops::Log::info() << "Grid config: " << *gridParams << std::endl;
    grid_ = atlas::Grid(*gridParams);
  } else if (iodaFile != boost::none) {
    oops::Log::info() << "Grid input file: " << *iodaFile << std::endl;

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
      if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
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

  // Default mask, set to 1 (true)
  gmask_ = functionSpace_.createField<int>(
    atlas::option::name("gmask") | atlas::option::levels(levels_));
  auto maskView = atlas::array::make_view<int, 2>(gmask_);
  for (atlas::idx_t jnode = 0; jnode < gmask_.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < gmask_.shape(1); ++jlevel) {
       maskView(jnode, jlevel) = 1;
    }
  }

  // Land-sea mask
  const boost::optional<bool> &landsea_mask = params.landsea_mask.value();
  if (landsea_mask != boost::none) {
    if (*landsea_mask) {
      // Lon/lat sizes
      size_t nlon = 0;
      size_t nlat = 0;

      // NetCDF IDs
      int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

      if (comm_.rank() == 0) {
        // Open NetCDF file
        if ((retval = nc_open("/home/benjamin/code/bundle/saber/test/testdata/landsea.nc", NC_NOWRITE, &ncid))) ERR(retval);

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
      }

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);

      // Broadcast coordinates and land-sea mask
      comm_.broadcast(lon.begin(), lon.end(), 0);
      comm_.broadcast(lat.begin(), lat.end(), 0);
      comm_.broadcast(lsm.begin(), lsm.end(), 0);
      double dlon = lon[1]-lon[0];
      double dlat = lat[1]-lat[0];

      // Lon/lat variables
      double zlon, zlat, dzlon, dzlat;
      int ilonm, ilonp, ilon, ilatm, ilatp, ilat;

      if (functionSpace_.type() == "StructuredColumns") {
        // StructuredColumns
        atlas::functionspace::StructuredColumns fs(functionSpace_);

        // Create local coordinates fieldset
        auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
        for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
          // Longitude
          zlon = lonlatView(jnode, 0);
          if (zlon < lon[0]) { 
            ilonm = nlon-1;
            ilonp = 0;
            dzlon = zlon-lon[ilonm]+360.0;
          } else if (zlon < lon[nlon-1]) {
            ilonm = (zlon-lon[0])/dlon;
            ilonp = ilonm+1;
            dzlon = zlon-lon[ilonm];
          } else {
            ilonm = nlon-1;
            ilonp = 0;
            dzlon = zlon-lon[ilonm];
          }
          if (dzlon < 0.5*dlon) {
            ilon = ilonm;
          } else {
            ilon = ilonp;
          }

          // Latitude
          zlat = lonlatView(jnode, 1);
          if (zlat < lat[0]) { 
            ilatm = 0;
            ilatp = 0;
            dzlat = 0.0;
          } else if (zlat < lat[nlat-1]) {
            ilatm = (zlat-lat[0])/dlat;
            ilatp = ilatm+1;
            dzlat = zlat-lat[ilatm];
          } else {
            ilatm = nlat-1;
            ilatp = nlat-1;
            dzlat = 0.0;
          }
          if (dzlat < 0.5*dlat) {
            ilat = ilatm;
          } else {
            ilat = ilatp;
          }

          // Get nearest neighbor value
          for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
            // Keep ocean points only
            if (lsm[ilat*nlon+ilon] == 0) {
              maskView(jnode, jlevel) = 1;
            } else {
              maskView(jnode, jlevel) = 0;
            }
          }
        }
      } else if (functionSpace_.type() == "NodeColumns") {
        // NodeColumns
        if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
// TODO(Benjamin): remove this line once ATLAS is upgraded to 0.29.0 everywhere
#if atlas_TRANS_FOUND
          // CubedSphere
          atlas::functionspace::CubedSphereNodeColumns fs(functionSpace_);

          // Create local coordinates fieldset
          auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
          for (atlas::idx_t jnode = 0; jnode < fs.lonlat().shape(0); ++jnode) {
            zlon = lonlatView(jnode, 0);
            zlat = lonlatView(jnode, 1);
          }
#else
          ABORT("TRANS required");
#endif
        } else {
          // Other NodeColumns
          atlas::functionspace::NodeColumns fs(functionSpace_);

          // Create local coordinates fieldset
          auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());
          for (atlas::idx_t jnode = 0; jnode < fs.lonlat().shape(0); ++jnode) {
            zlon = lonlatView(jnode, 0);
            zlat = lonlatView(jnode, 1);
          }
        }
      } else {
        ABORT(functionSpace_.type() + " function space not supported yet");
      }
    }
  }

  // Fill extra geometry fields
  extraFields_ = atlas::FieldSet();

  // Vertical unit
  atlas::Field vunit = functionSpace_.createField<double>(
    atlas::option::name("vunit") | atlas::option::levels(levels_));
  auto view = atlas::array::make_view<double, 2>(vunit);
  for (atlas::idx_t jnode = 0; jnode < vunit.shape(0); ++jnode) {
    for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
       view(jnode, jlevel) = vunit_[jlevel];
    }
  }
  extraFields_->add(vunit);

  // Geographical mask
  extraFields_->add(gmask_);

  // Halo mask
  if (grid_.name().compare(0, 1, std::string{"L"}) == 0) {
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
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
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
  atlas::Field gmask = (*other.extraFields())["gmask"];
  extraFields_->add(gmask);
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
  os << "Mask size: " << "TODO" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace quench
