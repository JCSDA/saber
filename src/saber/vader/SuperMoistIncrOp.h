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

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

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
///        specific_humidity, mass_content_of_cloud_ice_in_atmosphere_layer
///        and mass_content_of_cloud_liquid_water_in_atmosphere_layer. To do this
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
                   const std::vector<size_t> &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const atlas::FieldSet &,
                   const atlas::FieldSet &);
  virtual ~SuperMoistIncrOp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  oops::Variables activeVars_;
  std::unique_ptr<AirTemperature> exnerThetaToTemp_;
  std::unique_ptr<MoistIncrOp> MIO_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
