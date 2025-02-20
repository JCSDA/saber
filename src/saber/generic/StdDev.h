/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
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

  // Scaling parameter
  oops::Parameter<double> scaleFactorParam{"stddev scale factor",
                                           "multiplicative factor applied to StdDev block",
                                           1.0, this,
                                          {oops::exclusiveMinConstraint(0.)}};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}

  oops::Parameter<eckit::LocalConfiguration> scaling{"standard deviations",
                                                     eckit::LocalConfiguration(), this};
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
         const oops::FieldSet3D &,
         const oops::FieldSet3D &);
  virtual ~StdDev() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override;

  void write() const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  Parameters_ params_;
  bool readFromAtlas_;
  bool readFromModel_;
  double scaleFactor_;
  eckit::LocalConfiguration readConf_;
  std::unique_ptr<oops::FieldSet3D> stdDevFset_;
  bool writeToAtlas_;
  bool writeToModel_;
  eckit::LocalConfiguration writeConf_;
  std::map<std::string, double> scaling_;

  // Interative mean
  std::unique_ptr<oops::FieldSet3D> iterativeMean_;
  // Interative variance
  std::unique_ptr<oops::FieldSet3D> iterativeVar_;
  // Interative counter
  size_t iterativeN_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
