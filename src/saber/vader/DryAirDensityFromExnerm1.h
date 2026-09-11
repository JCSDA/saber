/*
 * (C) Crown Copyright 2024-2025 Met Office
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

#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class DryAirDensityFromExnerm1Parameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensityFromExnerm1Parameters, BlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "air_potential_temperature",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "dimensionless_exner_function_levels_minus_one",
    "dry_air_density_levels_minus_one",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    }});}

  const oops::Variables mandatoryStateVars() const override {return oops::Variables({
    "air_potential_temperature",
    "air_pressure_levels",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "dry_air_density_levels_minus_one",
    "height_above_mean_sea_level_levels",
    "height_above_mean_sea_level",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    });}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["dry_air_density_levels_minus_one"].getLevels();
    oops::Variables vars;
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    vars.push_back({"air_potential_temperature", conf});
    vars.push_back({"dimensionless_exner_function_levels_minus_one", conf});
    vars.push_back({"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["dry_air_density_levels_minus_one"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class DryAirDensityFromExnerm1 : public OuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::DryAirDensityFromExnerm1";}

  typedef DryAirDensityFromExnerm1Parameters Parameters_;

  DryAirDensityFromExnerm1(const oops::GeometryData &,
                           const oops::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);
  virtual ~DryAirDensityFromExnerm1();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;
  void directCalibration(const oops::FieldSets &) override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
