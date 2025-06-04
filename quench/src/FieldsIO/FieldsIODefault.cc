/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/FieldsIO/FieldsIODefault.h"

#include <vector>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "src/Geometry.h"

namespace quench {

// -----------------------------------------------------------------------------

static FieldsIOMaker<FieldsIODefault> makerDefault_("default");

// -----------------------------------------------------------------------------

void FieldsIODefault::read(const Geometry & geom,
                           const oops::Variables & vars,
                           const eckit::Configuration & conf,
                           atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Create variableSizes
  std::vector<size_t> variableSizes;
  for (const auto & var : vars) {
    variableSizes.push_back(var.getLevels());
  }

  // Update configuration
  eckit::LocalConfiguration updatedConf(conf);
  if (!updatedConf.has("latitude south to north")) {
    updatedConf.set("latitude south to north", geom.latSouthToNorth());
  }

  // Read fieldset
  util::readFieldSet(geom.getComm(),
                     geom.functionSpace(),
                     variableSizes,
                     vars.variables(),
                     updatedConf,
                     fset);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsIODefault::write(const Geometry & geom,
                            const eckit::Configuration & conf,
                            const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Update configuration
  eckit::LocalConfiguration updatedConf(conf);
  if (!updatedConf.has("latitude south to north")) {
    updatedConf.set("latitude south to north", geom.latSouthToNorth());
  }

  // Write fieldset
  util::writeFieldSet(geom.getComm(), updatedConf, fset);

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
