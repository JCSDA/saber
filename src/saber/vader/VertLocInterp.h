/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------
class VertLocInterpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VertLocInterpParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<atlas::idx_t> innerVerticalLevels{"inner vertical levels",
    "inner number of vertical levels", this};
  oops::Parameter<bool> reproduceBugStaggerDefn{
    "reproduce bug stagger definition",
    false, this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class VertLocInterp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::VertLocInterp";}

  typedef VertLocInterpParameters Parameters_;

  VertLocInterp(const oops::GeometryData &,
                const oops::Variables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const oops::FieldSet3D &,
                const oops::FieldSet3D &);
  virtual ~VertLocInterp() = default;

  const oops::GeometryData & innerGeometryData()
    const override {return outerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  Parameters_ params_;
  const oops::GeometryData & outerGeometryData_;
  oops::Variables outerVars_;
  oops::Variables activeVars_;
  oops::Variables innerVars_;
};

}  // namespace vader
}  // namespace saber
