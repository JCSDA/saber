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
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

namespace saber {
namespace interpolation {

class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  /// Function space
  oops::RequiredParameter<std::string> functionSpace{"function space", this};

  /// Grid
  oops::OptionalParameter<eckit::LocalConfiguration> grid{"grid", this};

  /// Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "equal_regions", this};

  /// Halo size
  oops::OptionalParameter<size_t> halo{"halo", this};
};

// -----------------------------------------------------------------------------
/// Geometry handles geometry for model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  typedef GeometryParameters Parameters_;

  static const std::string classname() {return "interpolation::Geometry";}

  Geometry(const Parameters_ &,
           const eckit::mpi::Comm & comm = eckit::mpi::comm());

  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  size_t halo_;
  atlas::Grid grid_;
  bool unstructuredGrid_;
  atlas::grid::Partitioner partitioner_;
  atlas::FunctionSpace functionSpace_;
};
// -----------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace saber
