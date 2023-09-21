/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/grid/Grid.h"

#include <memory>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/gsi/grid/Grid.interface.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

Grid::Grid(const eckit::mpi::Comm & comm, const eckit::Configuration & conf)
{
  oops::Log::trace() << classname() << "::Grid starting" << std::endl;
  util::Timer timer(classname(), "Grid");

  // Create grid
  gsi_grid_create_f90(keySelf_, conf, comm);

  // Get number of levels
  gsi_grid_get_levels_f90(keySelf_, gsiLevels_);

  // Create a functionspace for the GSI grid
  atlas::FieldSet gsiGridFieldSet = atlas::FieldSet();
  gsi_grid_set_atlas_lonlat_f90(keySelf_, gsiGridFieldSet.get());
  atlas::Field lonlat = gsiGridFieldSet["lonlat"];
  gsiGridFuncSpace_ = atlas::functionspace::PointCloud(lonlat);

  oops::Log::trace() << classname() << "::Grid done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Grid::~Grid() {
  oops::Log::trace() << classname() << "::~Grid starting" << std::endl;
  util::Timer timer(classname(), "~Grid");
  gsi_grid_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Grid done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Grid::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  gsi_grid_print_f90(keySelf_);
  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
