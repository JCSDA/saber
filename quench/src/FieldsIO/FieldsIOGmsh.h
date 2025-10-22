/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "src/FieldsIO/FieldsIOBase.h"

namespace quench {
  class Fields;

// -----------------------------------------------------------------------------
/// FieldsIOGmsh class

class FieldsIOGmsh : public FieldsIOBase {
 public:
  static const std::string classname()
    {return "quench::FieldsIOGmsh";}

  // Constructor/destructor
  explicit FieldsIOGmsh(const std::string & ioFormat)
    : FieldsIOBase(ioFormat) {}
  ~FieldsIOGmsh() = default;

  // Write
  void write(const eckit::Configuration &,
             const Fields &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace quench

