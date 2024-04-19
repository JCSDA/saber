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

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/bump/BUMP.h"
#include "saber/bump/BUMPParameters.h"

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

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  oops::Variables activeVars_;
  BUMPParameters bumpParams_;
  std::unique_ptr<BUMP> bump_;
  size_t memberIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
