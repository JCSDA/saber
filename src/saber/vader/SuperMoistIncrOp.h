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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/AirTemperature.h"
#include "saber/vader/MoistIncrOp.h"
#include "saber/vader/MoistIncrOpParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {
  class AirTemperature;
  class MoistIncrOp;

// -----------------------------------------------------------------------------
/// \brief a super saber block that splits total water ("qt") into
///        water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water (AKA specific_humidity),
///        cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water
///        and cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water. To do this
///        separation, the MoistIncrOp saber block requires the air_temperature increment.
///        The calibration also requires the air_temperature increment as input.
///        To make sure that the order of saber blocks is unaffected, this super
///        saber block has been created that uses the MoistIncrOp and AirTemperature
///        saber blocks.
class SuperMoistIncrOp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::SuperMoistIncrOp";}

  typedef SuperMoistIncrOpParameters Parameters_;

  SuperMoistIncrOp(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~SuperMoistIncrOp();

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
  const oops::Variables intermediateTempVars_;
  std::unique_ptr<MoistIncrOp> MIO_;
  std::unique_ptr<AirTemperature> exnerThetaToTemp_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
