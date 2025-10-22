/*
 * (C) Copyright 2024-     UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/VariableChange.h"

#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "src/Geometry.h"
#include "src/State.h"

namespace quench {

// -----------------------------------------------------------------------------

VariableChange::VariableChange(const eckit::Configuration & config,
                               const Geometry & geom)
  : geom_(geom) {
  oops::Log::trace() << classname() << "::VariableChange starting" << std::endl;
  oops::Log::trace() << classname() << "::VariableChange done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVar(State & x,
                               const oops::Variables & vars_out) const {
  oops::Log::trace() << classname() << "::changeVar starting" << std::endl;

  // Create FieldSet
  atlas::FieldSet fset;

  // State to FieldSet
  x.toFieldSet(fset);

  // Create any fields that do not exist
  for (auto & var : vars_out) {
    if (!fset.has(var.name())) {
      fset.add(geom_.functionSpace().createField<double>(
        atlas::option::name(var.name()) |
        atlas::option::levels(geom_.levels(var.name()))));
    }
  }

  // At this point all the fields are available but no actual variable transform
  // has been completed. If the use case arises they should be added here.
  // For example if x contains winds and vars_out contain stream function and
  // velocity potential, they should be computed. Fields in x that are not
  // in vars_out should be removed but this might require changes to saber
  // functionality

  // FieldSet to State
  x.fromFieldSet(fset);

  oops::Log::trace() << classname() << "::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x,
                                      const oops::Variables & vars_out) const {
  oops::Log::trace() << classname() << "::changeVarInverse starting" << std::endl;

  // Not implemented yet
  throw eckit::NotImplemented("changeVarInverse not implemented", Here());

  oops::Log::trace() << classname() << "::changeVarInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
