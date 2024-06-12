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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/MoistIncrOpParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class MoistIncrOp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::MoistIncrOp";}

  typedef MoistIncrOpParameters Parameters_;

  MoistIncrOp(const oops::GeometryData &,
              const oops::Variables &,
              const eckit::Configuration &,
              const Parameters_ &,
              const oops::FieldSet3D &,
              const oops::FieldSet3D &);
  virtual ~MoistIncrOp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void directCalibration(const oops::FieldSets &) override;


 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  atlas::FieldSet augmentedStateFieldSet_;
  const Parameters_ blockparams_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
