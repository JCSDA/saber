/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "src/LinearVariableChangeParameters.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Geometry;
  class State;
  class Increment;

// -----------------------------------------------------------------------------
/// quench linear change of variable

class LinearVariableChange: public util::Printable {
 public:
  typedef LinearVariableChangeParameters Parameters_;
  static const std::string classname() {return "quench::LinearVariableChange";}

  LinearVariableChange(const Geometry &, const Parameters_ &) {}

/// Perform linear transforms
  void changeVarTL(Increment &, const oops::Variables &) const {}
  void changeVarInverseTL(Increment &, const oops::Variables &) const {}
  void changeVarAD(Increment &, const oops::Variables &) const {}
  void changeVarInverseAD(Increment &, const oops::Variables &) const {}

  void changeVarTraj(const State &, const oops::Variables &) {}

 private:
  void print(std::ostream & os) const override {os << "LinearVariableChange";};
};
// -----------------------------------------------------------------------------

}  // namespace quench
