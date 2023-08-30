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
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/bump/BUMPParameters.h"
#include "saber/bump/lib/BUMP.h"


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
        const oops::FieldSet3D &,
        const size_t &);
  virtual ~NICAS();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead() override;

  void read() override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const atlas::FieldSet &) override;
  void iterativeCalibrationFinal() override;

  void dualResolutionSetup(const oops::GeometryData &) override;

  void write() const override;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const override;

 private:
  void print(std::ostream &) const override;
  BUMPParameters bumpParams_;
  oops::Variables activeVars_;
  std::unique_ptr<bump_lib::BUMP> bump_;
  size_t memberIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
