/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

namespace saber {

// -----------------------------------------------------------------------------
/// Base class for SABER block chains that have a self-adjoint central block
/// (ensemble and non-ensemble).
class SaberBlockChainBase {
 public:
  SaberBlockChainBase() = default;
  virtual ~SaberBlockChainBase() = default;

  virtual void randomize(atlas::FieldSet &) const = 0;
  virtual void multiply(atlas::FieldSet &) const = 0;
};

}  // namespace saber
