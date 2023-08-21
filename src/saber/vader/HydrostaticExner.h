/*
 * (C) Crown Copyright 2022 Met Office
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
class HydrostaticExnerParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<GpToHpCovarianceParameters>
    hydrostaticexnerParams{"covariance data", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "air_pressure_levels",
    "exner_levels_minus_one",
    "geostrophic_pressure_levels_minus_one",
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels",
    "unbalanced_pressure_levels_minus_one"});}
};

// -----------------------------------------------------------------------------
/// \brief  This saber block is here to do 3 jobs:
///         1) the vertical regression on geostrophic pressure
///         2) summing the result with unbalanced pressure to
///            create hydrostatic_pressure
///         3) converting hydrostatic pressure to exner pressure.
// TO DO: Marek - remove saber block when results identical to 3 new saber blocks
class HydrostaticExner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydrostaticExner";}

  typedef HydrostaticExnerParameters Parameters_;

  HydrostaticExner(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const atlas::FieldSet &,
                   const atlas::FieldSet &,
                   const util::DateTime &);
  virtual ~HydrostaticExner();

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
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
