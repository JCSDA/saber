/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/bump/BUMPParameters.h"
#include "saber/bump/lib/BUMP.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class VerticalBalanceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VerticalBalanceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<BUMPParameters> readParams{"read", this};
  oops::OptionalParameter<BUMPParameters> calibrationParams{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class VerticalBalance : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::bump::VerticalBalance";}

  typedef VerticalBalanceParameters Parameters_;

  VerticalBalance(const oops::GeometryData &,
                  const oops::Variables &,
                  const eckit::Configuration &,
                  const Parameters_ &,
                  const oops::FieldSet3D &,
                  const oops::FieldSet3D &);
  virtual ~VerticalBalance();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

  void read() override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const atlas::FieldSet &) override;
  void iterativeCalibrationFinal() override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  BUMPParameters bumpParams_;
  oops::Variables activeVars_;
  std::unique_ptr<bump_lib::BUMP> bump_;
  size_t memberIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
