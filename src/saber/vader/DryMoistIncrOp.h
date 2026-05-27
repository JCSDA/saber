/*
 * (C) Crown Copyright 2026 Met Office
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

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

class DryMoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DryMoistIncrOpParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "water_vapor_mixing_ratio_wrt_dry_air",
    "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_liquid_water_mixing_ratio_wrt_dry_air",
    "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
    "cloud_ice_mixing_ratio_wrt_dry_air",
    "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"
    }});}

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({
      "water_vapor_mixing_ratio_wrt_dry_air",
      "cloud_liquid_water_mixing_ratio_wrt_dry_air",
      "cloud_ice_mixing_ratio_wrt_dry_air",
      "total_water_mixing_ratio_wrt_dry_air"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["water_vapor_mixing_ratio_wrt_dry_air"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water", conf});
    vars.push_back({"cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water", conf});

    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["water_vapor_mixing_ratio_wrt_dry_air"],
                          outerVars["cloud_liquid_water_mixing_ratio_wrt_dry_air"],
                          outerVars["cloud_ice_mixing_ratio_wrt_dry_air"]
                         });
    return vars;
  }
};

// -----------------------------------------------------------------------------
/// \brief This saber block is here to copy hydrostatic_exner_levels into
///        dimensionless_exner_function_levels_minus_one and hydrostatic_pressure_levels into
///        air_pressure_levels

class DryMoistIncrOp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::DryMoistIncrOp";}

  typedef DryMoistIncrOpParameters Parameters_;

  DryMoistIncrOp(const oops::GeometryData &,
                 const oops::Variables &,
                 const eckit::Configuration &,
                 const Parameters_ &,
                 const oops::FieldSet3D &,
                 const oops::FieldSet3D &);
  virtual ~DryMoistIncrOp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D & fset) const override
    {inverseMultiply(fset);}
  void rightInverseMultiply(oops::FieldSet3D & fset) const override
    {inverseMultiply(fset);}
  void directCalibration(const oops::FieldSets &) override;

 private:
  void print(std::ostream &) const override;
  void inverseMultiply(oops::FieldSet3D &) const;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
