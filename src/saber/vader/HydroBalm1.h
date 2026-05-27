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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HydroBalm1Parameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydroBalm1Parameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({
        std::vector<std::string>{
        "dimensionless_exner_function_levels_minus_one",
        "hydrostatic_exner_levels",
        "virtual_potential_temperature"}});
  }

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({
        "air_pressure_levels",
        "hydrostatic_pressure_levels",
        "virtual_potential_temperature",
        "height_above_mean_sea_level_levels"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels =
      outerVars["dimensionless_exner_function_levels_minus_one"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    eckit::LocalConfiguration conf2;
    conf2.set("levels", modelLevels+1);
    oops::Variables vars;
    vars.push_back({"air_pressure_levels_minus_one", conf});
    vars.push_back({"hydrostatic_exner_levels", conf2});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["virtual_potential_temperature"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------

class HydroBalm1 : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydroBalm1";}

  typedef HydroBalm1Parameters Parameters_;

  HydroBalm1(const oops::GeometryData &,
           const oops::Variables &,
           const eckit::Configuration &,
           const Parameters_ &,
           const oops::FieldSet3D &,
           const oops::FieldSet3D &);
  virtual ~HydroBalm1();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;
  void rightInverseMultiply(oops::FieldSet3D &) const override;
  void directCalibration(const oops::FieldSets &) override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
