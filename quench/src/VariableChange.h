/*
 * (C) Copyright 2017-2021 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
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
  static const std::string classname()
    {return "quench::VariableChange";}

  // Constructor/destructor
  VariableChange(const eckit::Configuration &, const Geometry &) {}
  ~VariableChange() {}

  // Variable changes: direct and inverse
  void changeVar(State &,
                 const oops::Variables &) const
    {throw eckit::NotImplemented(Here());}
  void changeVarInverse(State &,
                        const oops::Variables &) const
    {throw eckit::NotImplemented(Here());}

 private:
  // Print
  void print(std::ostream & os) const override
    {os << "VariableChange";};
};

// -----------------------------------------------------------------------------

}  // namespace quench
