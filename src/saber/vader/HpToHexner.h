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

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HpToHexnerParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HpToHexnerParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "hydrostatic_exner_levels",
    "hydrostatic_pressure_levels"});}
};

// -----------------------------------------------------------------------------
/// \brief This saber block converts
///        hydrostatic pressure to hydrostatic exner pressure.

class HpToHexner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HpToHexner";}

  typedef HpToHexnerParameters Parameters_;

  HpToHexner(const oops::GeometryData &,
             const std::vector<size_t> &,
             const oops::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const atlas::FieldSet &,
             const atlas::FieldSet &,
             const util::DateTime &);
  virtual ~HpToHexner();

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