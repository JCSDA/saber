/*
 * (C) Crown Copyright 2025 Met Office
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

#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

class Hpm1ToHexnerExnerm1Parameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(Hpm1ToHexnerExnerm1Parameters, BlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "dimensionless_exner_function_levels_minus_one",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels_minus_one"}});}

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({"air_pressure_levels",
                            "dimensionless_exner_function_levels_minus_one",
                            "hydrostatic_exner_levels",
                            "hydrostatic_pressure_levels"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["dimensionless_exner_function_levels_minus_one"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"hydrostatic_pressure_levels_minus_one", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["hydrostatic_exner_levels"],
                          outerVars["dimensionless_exner_function_levels_minus_one"]});
    return vars;
  }
};

// -----------------------------------------------------------------------------
/// \brief This saber block is here to copy hydrostatic_exner_levels into
///        dimensionless_exner_function_levels_minus_one and hydrostatic_pressure_levels into
///        air_pressure_levels

class Hpm1ToHexnerExnerm1 : public OuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::Hpm1ToHexnerExnerm1";}

  typedef Hpm1ToHexnerExnerm1Parameters Parameters_;

  Hpm1ToHexnerExnerm1(const oops::GeometryData &,
                     const oops::Variables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const oops::FieldSet3D &,
                     const oops::FieldSet3D &);
  virtual ~Hpm1ToHexnerExnerm1();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;
  void rightInverseMultiply(oops::FieldSet3D &) const override;

  void directCalibration(const oops::FieldSets & fset) override;

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
