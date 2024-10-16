/*
 * (C) Crown Copyright 2022-2024 Met Office
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

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class DryAirDensityParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryAirDensityParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "dry_air_density_levels_minus_one",
    "air_pressure_levels",
    "air_potential_temperature",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"}});}

  const oops::Variables mandatoryStateVars() const override {return oops::Variables({
    "dry_air_density_levels_minus_one",
    "air_pressure_levels_minus_one",
    "air_potential_temperature",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["dry_air_density_levels_minus_one"].getLevels();
    oops::Variables vars;
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels + 1);
    vars.push_back({"air_pressure_levels", conf});
    conf.set("levels", modelLevels);
    vars.push_back({"air_potential_temperature", conf});
    vars.push_back({"water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["dry_air_density_levels_minus_one"]});
    return vars;
  }

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
    if (outerVars.has("air_pressure_levels_minus_one")) {
      throw eckit::UserError("air_pressure_levels_minus_one is a "
                             "temporary variable of mo_dry_air_density "
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels = outerVars["dry_air_density_levels_minus_one"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables tempVars;
    tempVars.push_back({"air_pressure_levels_minus_one", conf});
    return tempVars;
  }
};

// -----------------------------------------------------------------------------

class DryAirDensity : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::DryAirDensity";}

  typedef DryAirDensityParameters Parameters_;

  DryAirDensity(const oops::GeometryData &,
                const oops::Variables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const oops::FieldSet3D &,
                const oops::FieldSet3D &);
  virtual ~DryAirDensity();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  const oops::Variables intermediateTempVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
