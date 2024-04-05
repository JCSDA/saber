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
    "dry_air_density_levels_minus_one",
    "air_pressure_levels",
    "potential_temperature",
    "specific_humidity",
    "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
    "mass_content_of_cloud_ice_in_atmosphere_layer"});}

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"air_pressure_levels",
                          "potential_temperature",
                          "specific_humidity",
                          "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
                          "mass_content_of_cloud_ice_in_atmosphere_layer"});
    const int modelLevels = outerVars.getLevels("dry_air_density_levels_minus_one");
    vars.addMetaData("air_pressure_levels", "levels", modelLevels + 1);
    vars.addMetaData("potential_temperature", "levels", modelLevels);
    vars.addMetaData("specific_humidity", "levels", modelLevels);
    vars.addMetaData("mass_content_of_cloud_liquid_water_in_atmosphere_layer", "levels",
                     modelLevels);
    vars.addMetaData("mass_content_of_cloud_ice_in_atmosphere_layer", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"dry_air_density_levels_minus_one"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }

  oops::Variables intermediateTempVars(const oops::Variables& outerVars) const {
    oops::Variables tempVars({"air_pressure_levels_minus_one"});
    if (outerVars.has("air_pressure_levels_minus_one")) {
      throw eckit::UserError("air_pressure_levels_minus_one is a "
                             "temporary variable of mo_dry_air_density "
                             " and should not be an outer variable of this block.",
                             Here());
    }
    const int modelLevels = outerVars.getLevels("dry_air_density_levels_minus_one");
    tempVars.addMetaData("air_pressure_levels_minus_one", "levels", modelLevels);
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
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
