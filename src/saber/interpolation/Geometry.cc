/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/Geometry.h"

#include <cmath>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

#include "saber/interpolation/GsiGrid.h"

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm)
  : comm_(comm), halo_(1)
{
  if (config.has(GsiGridKey) != config.has(GsiPartitionerKey)) {
    throw eckit::BadParameter("Must specify GSI-matching grid AND partitioner, OR neither");
  }

  if (config.has(GsiGridKey)) {
    // Here we've requested a particular grid and MPI partitioner setup. This can be used to
    // recreate the GSI grid and MPI partition, though not the GSI grid point ordering on each MPI
    // task. This is useful for the GSI saber block, but could be moved into the generic OOPS
    // utilities if desired elsewhere. We do this with custom code because GSI orders its grid
    // points (and MPI partitions) from south-to-north whereas atlas orders its grid points
    // (and partitions) from north-to-south.
    setupGsiMatchingGrid(config, comm_, grid_, functionSpace_, fieldSet_);
  } else {
    // Generic case
    atlas::Mesh mesh;
    util::setupFunctionSpace(comm_, config, grid_, partitioner_, mesh, functionSpace_, fieldSet_);
  }

  if (config.has("halo")) {
    halo_ = config.getUnsigned("halo");
  }

  // Print summary
  this->print(oops::Log::info());
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
