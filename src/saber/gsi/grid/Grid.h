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

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/gsi/grid/Grid.interface.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class GridParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GridParameters, Parameters)
};

// -------------------------------------------------------------------------------------------------

class Grid {
 public:
  static const std::string classname() {return "saber::gsi::Grid";}

  // Constructor & destructor
  Grid(const eckit::mpi::Comm &, const eckit::Configuration &);
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

}  // namespace gsi
}  // namespace saber
