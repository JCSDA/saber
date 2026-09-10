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

#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HpToHexnerParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HpToHexnerParameters, BlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    std::vector<std::string>{
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels"}});}

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables(std::vector<std::string>{"hydrostatic_exner_levels",
                            "hydrostatic_pressure_levels"});
  }
};

// -----------------------------------------------------------------------------
/// \brief This saber block converts
///        hydrostatic pressure to hydrostatic exner pressure.

class HpToHexner : public OuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HpToHexner";}

  typedef HpToHexnerParameters Parameters_;

  HpToHexner(const oops::GeometryData &,
             const oops::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const oops::FieldSet3D &,
             const oops::FieldSet3D &);
  virtual ~HpToHexner();

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
  oops::Variables activeVars_;
  oops::FieldSet3D xb_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
