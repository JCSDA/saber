/*
 * (C) Copyright 2024 Meteorlogisk Institutt
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
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class FakeLevelsParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FakeLevelsParametersBase, oops::Parameters)

 public:
  // Number of fake levels
  oops::OptionalParameter<size_t> nz{"number of fake levels", this};

  // Lowest fake level
  oops::OptionalParameter<double> lowestFakeLevel{"lowest fake level", this};

  // Highest fake level
  oops::OptionalParameter<double> highestFakeLevel{"highest fake level", this};

  // Explicit fake levels
  oops::OptionalParameter<std::vector<double>> fakeLevels{"fake levels", this};

  // Input model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputModelFilesConf{
    "input model files", this};

  // Output model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> outputModelFilesConf{
    "output model files", this};

  // Scalar vertical support
  oops::OptionalParameter<double> rvFromYaml{"vertical length-scale", this};
};

// -----------------------------------------------------------------------------

class FakeLevelsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(FakeLevelsParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<FakeLevelsParametersBase> read{"read", this};
  oops::OptionalParameter<FakeLevelsParametersBase> calibration{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class FakeLevels : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::generic::FakeLevels";}

  typedef FakeLevelsParameters     Parameters_;
  typedef FakeLevelsParametersBase ParametersBase_;

  FakeLevels(const oops::GeometryData &,
             const oops::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const oops::FieldSet3D &,
             const oops::FieldSet3D &);
  virtual ~FakeLevels() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return gdata_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
    override;

 private:
  const util::DateTime validTime_;
  const oops::GeometryData & gdata_;
  const eckit::mpi::Comm & comm_;
  oops::Variables outerVars_;
  oops::Variables activeVars_;
  const std::string suffix_;
  ParametersBase_ params_;
  eckit::LocalConfiguration fieldsMetaData_;
  size_t nz_;
  oops::Variables innerVars_;
  std::unique_ptr<oops::FieldSet3D> rv_;
  std::unique_ptr<oops::FieldSet3D> weight_;

  // Utilities
  eckit::LocalConfiguration getFileConf(const eckit::mpi::Comm &,
                                        const eckit::Configuration &) const;

  void print(std::ostream &) const override;
};

}  // namespace generic
}  // namespace saber
