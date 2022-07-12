/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_LINEARVARCHANGE_H_
#define QUENCH_LINEARVARCHANGE_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "quench/LinearVarChangeParams.h"

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

class LinearVarChange: public util::Printable {
 public:
  typedef LinearVarChangeParams Parameters_;
  static const std::string classname() {return "quench::LinearVarChange";}

  LinearVarChange(const Geometry &, const Parameters_ &) {}

/// Perform linear transforms
  void changeVarTL(Increment &, const oops::Variables &) const {}
  void changeVarInverseTL(Increment &, const oops::Variables &) const {}
  void changeVarAD(Increment &, const oops::Variables &) const {}
  void changeVarInverseAD(Increment &, const oops::Variables &) const {}

  void changeVarTraj(const State &, const oops::Variables &) {}

 private:
  void print(std::ostream &) const override {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
#endif  // QUENCH_LINEARVARCHANGE_H_
