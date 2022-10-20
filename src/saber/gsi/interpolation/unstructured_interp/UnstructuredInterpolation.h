/*
 * (C) Copyright 2020- UCAR
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

#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace saber {
namespace gsi {

// -----------------------------------------------------------------------------

/*! \brief Interface for Unstructured interpolation
 *
 */

class UnstructuredInterpolation : private util::ObjectCounter<UnstructuredInterpolation> {
 public:
  static const std::string classname() {return "oops::UnstructuredInterpolation";}

  UnstructuredInterpolation(const eckit::Configuration &,
                            const atlas::FunctionSpace &,
                            const atlas::FunctionSpace &,
                            const std::vector<std::string> &,
                            const eckit::mpi::Comm & = oops::mpi::world());
  ~UnstructuredInterpolation();

  void apply(const atlas::Field &, atlas::Field &);
  void apply(const atlas::FieldSet &, atlas::FieldSet &);

  void apply_ad(const atlas::Field &, atlas::Field &);
  void apply_ad(const atlas::FieldSet &, atlas::FieldSet &);

  int write(const eckit::Configuration &);

 private:
  int keyUnstructuredInterpolator_;
  const atlas::FunctionSpace *in_fspace_;
  const atlas::FunctionSpace *out_fspace_;
  std::vector<std::string> activeVars_;
  void print(std::ostream &) const;
};

// -----------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
