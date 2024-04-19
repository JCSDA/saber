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
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/bump/BUMP.h"
#include "saber/bump/BUMPParameters.h"


namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class NICASParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(NICASParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<BUMPParameters> readParams{"read", this};
  oops::OptionalParameter<BUMPParameters> calibrationParams{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class NICAS : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::bump::NICAS";}

  typedef NICASParameters Parameters_;

  NICAS(const oops::GeometryData &,
        const oops::Variables &,
        const eckit::Configuration &,
        const Parameters_ &,
        const oops::FieldSet3D &,
        const oops::FieldSet3D &);
  virtual ~NICAS();

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  void dualResolutionSetup(const oops::GeometryData &) override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override;

  size_t ctlVecSize() const override {return bump_->getCvSize();}
  void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const override;

 private:
  void print(std::ostream &) const override;
  oops::Variables activeVars_;
  BUMPParameters bumpParams_;
  std::unique_ptr<BUMP> bump_;
  size_t memberIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
