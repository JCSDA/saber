/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class VerticalBalanceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMPParameters> bumpParams{"bump", this};
  oops::Variables mandatoryActiveVars() const {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class VerticalBalance : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::bump::VerticalBalance";}

  typedef VerticalBalanceParameters Parameters_;

  VerticalBalance(const oops::GeometryData &,
                  const std::vector<size_t> &,
                  const oops::Variables &,
                  const Parameters_ &,
                  const atlas::FieldSet &,
                  const atlas::FieldSet &,
                  const std::vector<atlas::FieldSet> &);
  virtual ~VerticalBalance();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  std::unique_ptr<BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
