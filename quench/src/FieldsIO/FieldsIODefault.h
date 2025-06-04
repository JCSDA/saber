/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"

#include "src/FieldsIO/FieldsIOBase.h"

namespace quench {
  class Geometry;

// -----------------------------------------------------------------------------

class FieldsIODefault : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quench::FieldsIODefault";}

  // Constructor/destructor
  explicit FieldsIODefault(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIODefault() = default;

  // Read
  void read(const Geometry &,
            const oops::Variables &,
            const eckit::Configuration &,
            atlas::FieldSet &) const override;

  // Write
  void write(const Geometry &,
             const eckit::Configuration &,
             const atlas::FieldSet &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace quench
