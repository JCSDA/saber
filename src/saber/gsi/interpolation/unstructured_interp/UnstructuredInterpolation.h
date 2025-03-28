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

#include "oops/base/Variables.h"
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

  UnstructuredInterpolation(const eckit::mpi::Comm &,
                            const eckit::Configuration &,
                            const atlas::FunctionSpace &,
                            const atlas::FunctionSpace &,
                            const oops::Variables &);
  ~UnstructuredInterpolation();

  void apply(const atlas::Field &, atlas::Field &);
  void apply(atlas::FieldSet &);

  void applyAD(const atlas::Field &, atlas::Field &);
  void applyAD(atlas::FieldSet &);

  int write(const eckit::Configuration &);

 private:
  int keyUnstructuredInterpolator_;
  const atlas::FunctionSpace innerFuncSpace_;
  const atlas::FunctionSpace outerFuncSpace_;
  oops::Variables activeVars_;
  void print(std::ostream &) const;
};

// -----------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
