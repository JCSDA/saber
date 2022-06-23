/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/abor1_cpp.h"

#include "saber/gsi/grid/GSI_Grid.interface.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class GridParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GridParameters, SaberBlockParametersBase)

 public:
  // File containing grid and coefficients
  oops::RequiredParameter<std::string> GSIFile{"gsi error covariance file", this};
  oops::RequiredParameter<std::string> GSINML{"gsi berror namelist file", this};

  // Handle vertical top-2-bottom and vice-verse wrt to GSI
  oops::Parameter<bool> vflip{"flip vertical grid", true, this};

  // Processor layout
  oops::Parameter<size_t> layoutx{"processor layout x direction", 1, this};
  oops::Parameter<size_t> layouty{"processor layout y direction", 1, this};

  // Debugging mode
  oops::Parameter<bool> debugMode{"debugging mode", false, this};
  oops::Parameter<bool> bypassGSI{"debugging bypass gsi", false, this};
  oops::Parameter<bool> bypassGSIbe{"debugging deep bypass gsi B error", false, this};
};

// -------------------------------------------------------------------------------------------------

class Grid {
 public:
  static const std::string classname() {return "saber::gsi::Grid";}

  typedef GridParameters Parameters_;

  // Constructor & destructor
  Grid(const eckit::mpi::Comm &, const Parameters_ &);
  ~Grid();

  // Accessor functions
  int levels() {return gsiLevels_;}
  const atlas::FunctionSpace & functionSpace() const {return gsiGridFuncSpace_;}
  atlas::FunctionSpace & functionSpace() {return gsiGridFuncSpace_;}

 private:
  void print(std::ostream &) const;
  // Fortran LinkedList key
  GridKey keySelf_;
  // Function spaces
  atlas::FunctionSpace gsiGridFuncSpace_;
  // Number of levels
  int gsiLevels_;
};

// -------------------------------------------------------------------------------------------------

Grid::Grid(const eckit::mpi::Comm & comm, const Parameters_ & params)
{
  oops::Log::trace() << classname() << "::Grid starting" << std::endl;
  util::Timer timer(classname(), "Grid");

  // Create grid
  gsi_grid_create_f90(keySelf_, params.toConfiguration(), comm);

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
