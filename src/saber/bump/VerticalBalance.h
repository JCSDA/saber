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
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class VerticalBalanceParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceParameters, SaberOuterBlockParametersBase)

 public:
  oops::RequiredParameter<BUMPParameters> bumpParams{"bump", this};
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

  const oops::GeometryData & inputGeometryData() const override {return inputGeometryData_;}
  const oops::Variables & inputVars() const override {return inputVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & inputGeometryData_;
  oops::Variables inputVars_;
  std::unique_ptr<BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
