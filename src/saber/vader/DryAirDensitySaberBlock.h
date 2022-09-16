/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
class DryAirDensitySaberBlockParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensitySaberBlockParameters, SaberOuterBlockParametersBase)
 public:
};

// -----------------------------------------------------------------------------
class DryAirDensitySaberBlock : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::DryAirDensitySaberBlock";}

  typedef DryAirDensitySaberBlockParameters Parameters_;

  DryAirDensitySaberBlock(const eckit::mpi::Comm &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const eckit::Configuration &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);
  virtual ~DryAirDensitySaberBlock();

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace saber
