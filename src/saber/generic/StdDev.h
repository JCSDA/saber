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
#include "oops/util/parameters/Parameters.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class StdDevReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StdDevReadParameters, Parameters)

 public:
  // ATLAS standard-deviation file
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFileConf{"atlas file", this};
  // Model standard-deviation file
  oops::OptionalParameter<eckit::LocalConfiguration> modelFileConf{"model file", this};
};

// -----------------------------------------------------------------------------

class StdDevWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StdDevWriteParameters, Parameters)

 public:
  // ATLAS standard-deviation file
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFileConf{"write to atlas file", this};
  // Model standard-deviation file
  oops::OptionalParameter<eckit::LocalConfiguration> modelFileConf{"write to model file", this};
};

// -----------------------------------------------------------------------------

class StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(StdDevParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<StdDevReadParameters> readParams{"read", this};

  // Calibration of block parameters
  oops::OptionalParameter<StdDevWriteParameters> calibrationParams{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class StdDev : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::StdDev";}

  typedef StdDevParameters Parameters_;

  StdDev(const oops::GeometryData &,
         const oops::Variables &,
         const eckit::Configuration &,
         const Parameters_ &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const util::DateTime &);
  virtual ~StdDev() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead() override;

  void read() override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const atlas::FieldSet &) override;
  void iterativeCalibrationFinal() override;

  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const override;

  void write() const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  Parameters_ params_;
  bool readFromAtlas_;
  bool readFromModel_;
  eckit::LocalConfiguration readConf_;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> inputs_;
  atlas::FieldSet stdDevFset_;
  bool writeToAtlas_;
  bool writeToModel_;
  eckit::LocalConfiguration writeConf_;

  // Interative mean
  atlas::FieldSet iterativeMean_;
  // Interative variance
  atlas::FieldSet iterativeVar_;
  // Interative counter
  size_t iterativeN_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
