/*
 * (C) Crown Copyright 2023 Met Office
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

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/PressureParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------
/// \brief This saber block is here is to create "hydrostatic_pressure"
///        from "geostrophic pressure" and "unbalanced pressure"
///        Vertical regression is applied to "geostrophic pressure"
///        as part of this calculation.

class GpToHp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::GpToHp";}

  typedef GpToHpParameters Parameters_;

  GpToHp(const oops::GeometryData &,
         const oops::Variables &,
         const eckit::Configuration &,
         const Parameters_ &,
         const oops::FieldSet3D &,
         const oops::FieldSet3D &);
  virtual ~GpToHp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;
  void read() override;
  void directCalibration(const oops::FieldSets & fset) override;
  void write() const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  Parameters_ params_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

}  // namespace vader
}  // namespace saber
