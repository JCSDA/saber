/*
 * (C) Copyright 2023  UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <unordered_map>

#include "atlas/field.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Geometry;
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// LinearVariableChange class

class LinearVariableChange: public util::Printable,
                            private util::ObjectCounter<LinearVariableChange> {
 public:
  static const std::string classname()
    {return "quench::LinearVariableChange";}

  // Constructor/destructor
  LinearVariableChange(const Geometry &,
                       const eckit::Configuration &);
  ~LinearVariableChange()
    {}

  // Linear variable changes: TL, inverseTL, AD and inverseAD
  void changeVarTL(Increment &,
                   const oops::Variables &) const;
  void changeVarInverseTL(Increment &,
                          const oops::Variables &) const;
  void changeVarAD(Increment &,
                   const oops::Variables &) const;
  void changeVarInverseAD(Increment &,
                          const oops::Variables &) const;

  // Trajectory setup
  void changeVarTraj(const State &,
                     const oops::Variables &)
    {}

 private:
  // Print
  void print(std::ostream & os) const override
    {os << "LinearVariableChange";}

  // Map from output to input variables
  std::unordered_map<std::string, std::string> map_;

  // Multiplicative factor
  atlas::FieldSet multiplierFset_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
