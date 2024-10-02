/*
 * (C) Copyright 2023 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Printable.h"

#include "src/Geometry.h"

namespace quench {

// -------------------------------------------------------------------------------------------------

class ModelData : public util::Printable {
 public:
  static const std::string classname()
    {return "quench::ModelData";}

  // Constructor/destructor
  explicit ModelData(const Geometry & geometry) : modelData_(geometry.modelData()) {}
  ~ModelData() {}

  // Model data accessor
  const eckit::LocalConfiguration modelData() const {return modelData_;}

 private:
  // Print
  void print(std::ostream & os) const
    {os << "quench::ModelData::modelData(): " << modelData();}

  // Model data
  eckit::LocalConfiguration modelData_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace quench
