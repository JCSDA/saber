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

/// This is a helper class for setting up the SABER interpolation blocks -- it looks similar to a
/// model Geometry class (in particular, quench's), but is not fully conforming to the oops
/// interface. Its purpose is to help convert yaml configs into the desired atlas data that must be
/// passed into the interpolator.
class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "interpolation::Geometry";}

  explicit Geometry(const eckit::Configuration &,
                    const eckit::mpi::Comm & comm = eckit::mpi::comm());

  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  const atlas::FieldSet & fields() const {return fieldSet_;}

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  atlas::Grid grid_;
  atlas::grid::Partitioner partitioner_;
  size_t halo_;
  atlas::FunctionSpace functionSpace_;
  atlas::FieldSet fieldSet_;
};

}  // namespace interpolation
}  // namespace saber
