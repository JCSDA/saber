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

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace saber {
namespace interpolation {
// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), halo_(1)
{
  atlas::Mesh mesh;
  util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh, functionSpace_, fieldSet_);

  if (config.has("halo")) {
    halo_ = config.getUnsigned("halo");
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
  if (partitioner_) {
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
