/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "src/Geometry.h"
#include "src/State.h"

namespace quench {
  class State;

// -----------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname() {return "quench::VariableChange";}

  VariableChange(const eckit::Configuration &, const Geometry &) {}
  ~VariableChange() {}

/// Perform transforms
  void changeVar(State &, const oops::Variables &) const
    {throw eckit::NotImplemented(Here());}
  void changeVarInverse(State &, const oops::Variables &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream & os) const override {os << "VariableChange";};
};

// -----------------------------------------------------------------------------

}  // namespace quench
