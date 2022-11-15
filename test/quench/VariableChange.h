/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "quench/Geometry.h"
#include "quench/State.h"
#include "quench/VariableChangeParameters.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace quench {
  class State;

// -----------------------------------------------------------------------------
/// quench change of variable

class VariableChange : public util::Printable {
 public:
  typedef VariableChangeParameters Parameters_;
  static const std::string classname() {return "quench::VariableChange";}

  VariableChange(const Parameters_ &, const Geometry &) {}

/// Perform transforms
  void changeVar(State &, const oops::Variables &) const {}
  void changeVarInverse(State &, const oops::Variables &) const {}

 private:
  void print(std::ostream & os) const override {os << "VariableChange";};
};
// -----------------------------------------------------------------------------

}  // namespace quench
