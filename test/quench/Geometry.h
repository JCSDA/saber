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
  oops::RequiredParameter<eckit::LocalConfiguration> grid{"grid", this};

  /// Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "checkerboard", this};

  /// Number of levels
  oops::Parameter<size_t> levels{"levels", 1, this};

  /// Vertical unit
  oops::OptionalParameter<std::vector<double>> vunit{"vunit", this};
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
  atlas::Grid * atlasGrid() const {return atlasGrid_.get();}
  atlas::FunctionSpace * atlasFunctionSpace() const
    {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  size_t levels() const {return levels_;}
  std::vector<double> vunit() const {return vunit_;}

  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const {}

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  eckit::LocalConfiguration gridConfig_;
  std::unique_ptr<atlas::Grid> atlasGrid_;
  atlas::Mesh atlasMesh_;
  std::unique_ptr<atlas::FunctionSpace> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  size_t levels_;
  std::vector<double> vunit_;
};
// -----------------------------------------------------------------------------

}  // namespace quench

#endif  // QUENCH_GEOMETRY_H_
