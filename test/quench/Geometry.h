/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_GEOMETRY_H_
#define QUENCH_GEOMETRY_H_

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

  /// Grid input file (NetCDF format)
  oops::OptionalParameter<std::string> gridInput{"grid input file", this};

  /// Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "trans", this};

  /// Number of levels
  oops::Parameter<size_t> levels{"levels", 1, this};

  /// Vertical unit
  oops::OptionalParameter<std::vector<double>> vunit{"vunit", this};

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
  std::vector<double> vunit() const {return vunit_;}
  size_t halo() const {return halo_;}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const {}

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  atlas::Grid grid_;
  std::string gridType_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::FunctionSpace functionSpace_;
  atlas::FieldSet extraFields_;
  size_t levels_;
  std::vector<double> vunit_;
  size_t halo_;
};
// -----------------------------------------------------------------------------

}  // namespace quench

#endif  // QUENCH_GEOMETRY_H_
