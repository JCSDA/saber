/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/Geometry.h"

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

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace saber {
namespace interpolation {
// -----------------------------------------------------------------------------

Geometry::Geometry(const Parameters_ & params,
                   const eckit::mpi::Comm & comm) :
  comm_(comm) {
  // Halo
  const auto & halo = params.halo.value();
  if (halo != boost::none) {
    halo_ = *halo;
  } else {
    halo_ = 1;
  }

  // Setup grid
  const auto & gridParams = params.grid.value();
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
  } else {
    throw eckit::UserError("Grid or grid input file required", Here());
  }

  if (!unstructuredGrid_) {
    // Setup partitioner
    partitioner_ = atlas::grid::Partitioner(params.partitioner.value());
  }

  if (params.functionSpace.value() == "StructuredColumns") {
    // StructuredColumns

    // Setup distribution
    atlas::grid::Distribution distribution = atlas::grid::Distribution(grid_, partitioner_);

    // Setup function space
    functionSpace_ = atlas::functionspace::StructuredColumns(grid_, distribution,
                     atlas::option::halo(halo_));

  } else if (params.functionSpace.value() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      const atlas::Mesh mesh = atlas::MeshGenerator("cubedsphere_dual").generate(grid_);
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(mesh);
    } else {
      if (comm_.size() == 1) {
        // NodeColumns
        const atlas::Mesh mesh = atlas::MeshGenerator("delaunay").generate(grid_);
        functionSpace_ = atlas::functionspace::NodeColumns(mesh);
      } else {
        throw eckit::FunctionalityNotSupported("NodeColumns function space on multiple PEs not"
          " supported yet", Here());
      }
    }
  } else if (params.functionSpace.value() == "PointCloud") {
    // Setup function space
    functionSpace_ = atlas::functionspace::PointCloud(grid_);
  } else {
    throw eckit::FunctionalityNotSupported(params.functionSpace.value() +
      " function space not supported yet", Here());
  }

  // Print summary
  this->print(oops::Log::info());
}
// -------------------------------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool includeHalo) const {
  const auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  const auto ghost = atlas::array::make_view<int, 1>(functionSpace_.ghost());

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
  os << prefix <<  "Interpolation geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (!unstructuredGrid_) {
    os << prefix << "Partitioner:" << std::endl;
    os << prefix << "- type: " << partitioner_.type() << std::endl;
  }
  os << prefix << "Function space:" << std::endl;
  os << prefix << "- type: " << functionSpace_.type() << std::endl;
  os << prefix << "- halo: " << halo_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace interpolation
}  // namespace saber
