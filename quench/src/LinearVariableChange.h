/*
 * (C) Copyright 2023  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "atlas/field.h"

#include "oops/util/Printable.h"

#include "src/Geometry.h"
#include "src/Increment.h"
#include "src/State.h"

namespace quench {

// -----------------------------------------------------------------------------
/// quench linear change of variable (simple pointwise multiplication with a
/// field loaded from file)

class LinearVariableChange: public util::Printable {
 public:
  static const std::string classname() {return "quench::LinearVariableChange";}

  LinearVariableChange(const Geometry &, const eckit::Configuration &);
  ~LinearVariableChange();

/// Perform linear transforms
  void changeVarTL(Increment &, const oops::Variables &) const;
  void changeVarInverseTL(Increment &, const oops::Variables &) const;
  void changeVarAD(Increment &, const oops::Variables &) const;
  void changeVarInverseAD(Increment &, const oops::Variables &) const;

  void changeVarTraj(const State &, const oops::Variables &) {}

 private:
  void print(std::ostream & os) const override {os << "LinearVariableChange";};
  atlas::FieldSet fset_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
