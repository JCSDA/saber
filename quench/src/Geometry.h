/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;
}

namespace quench {

class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  /// Function space
  oops::RequiredParameter<std::string> functionSpace{"function space", this};

  /// Grid
  oops::OptionalParameter<eckit::LocalConfiguration> grid{"grid", this};

  /// IODA input file (NetCDF format) to get longitude/latitude
  oops::OptionalParameter<std::string> iodaFile{"ioda file", this};

  /// Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "equal_regions", this};

  /// Number of levels
  oops::Parameter<size_t> levels{"levels", 1, this};

  /// Corresponding level for 2D variables (first or last)
  oops::Parameter<std::string> lev2d{"lev2d", "first", this};

  /// Vertical unit
  oops::OptionalParameter<std::vector<double>> vunit{"vunit", this};

  /// Mask type
  oops::Parameter<std::string> mask_type{"mask type", "none", this};

  /// Mask path
  oops::Parameter<std::string> mask_path{"mask path", "../quench/data/landsea.nc", this};

  /// Halo size
  oops::OptionalParameter<size_t> halo{"halo", this};
};

// -----------------------------------------------------------------------------
/// Geometry handles geometry for  model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  typedef GeometryParameters Parameters_;

  static const std::string classname() {return "quench::Geometry";}

  Geometry(const Parameters_ &,
           const eckit::mpi::Comm & comm = oops::mpi::world());
  Geometry(const Geometry &);

  const eckit::mpi::Comm & getComm() const {return comm_;}
  const atlas::Grid grid() const {return grid_;}
  const std::string gridType() const {return gridType_;}
  const atlas::grid::Partitioner partitioner() const {return partitioner_;}
  const atlas::Mesh mesh() const {return mesh_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  atlas::FunctionSpace & functionSpace() {return functionSpace_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}
  size_t levels() const {return levels_;}

  size_t variableSize(const std::string &) const;
  size_t maskLevel(const std::string &, const size_t &) const;
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const {}
  bool levelsAreTopDown() const {return true;}

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  size_t halo_;
  atlas::Grid grid_;
  std::string gridType_;
  bool unstructuredGrid_;
  atlas::grid::Partitioner partitioner_;
  atlas::grid::Distribution distribution_;
  atlas::Mesh mesh_;
  atlas::Field gmask_;
  double gmaskSize_;
  atlas::FunctionSpace functionSpace_;
  atlas::FieldSet extraFields_;
  size_t levels_;
  std::string lev2d_;
  std::vector<double> vunit_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
