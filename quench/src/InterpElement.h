/*
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace quench {

// -----------------------------------------------------------------------------

class InterpElement {
 public:
  // Constructor
  explicit InterpElement(const std::vector<std::pair<size_t, double>> & operations) :
    operations_(operations) {}

  // Accessors
  const std::vector<std::pair<size_t, double>> & operations() const {return operations_;}

 private:
  // Operations
  std::vector<std::pair<size_t, double>> operations_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
