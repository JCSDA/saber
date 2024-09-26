/*
 * (C) Copyright 2024 Meteorologisk Institutt
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

class ShadowLevelsParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ShadowLevelsParametersBase, oops::Parameters)

 public:
  // Number of shadow levels
  oops::OptionalParameter<size_t> nz{"number of shadow levels", this};

  // Lowest shadow level
  oops::OptionalParameter<double> lowestShadowLevel{"lowest shadow level", this};

  // Highest shadow level
  oops::OptionalParameter<double> highestShadowLevel{"highest shadow level", this};

  // Explicit shadow levels
  oops::OptionalParameter<std::vector<double>> shadowLevels{"shadow levels", this};

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

class ShadowLevelsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ShadowLevelsParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<ShadowLevelsParametersBase> read{"read", this};
  oops::OptionalParameter<ShadowLevelsParametersBase> calibration{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class ShadowLevels : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::generic::ShadowLevels";}

  typedef ShadowLevelsParameters     Parameters_;
  typedef ShadowLevelsParametersBase ParametersBase_;

  ShadowLevels(const oops::GeometryData &,
               const oops::Variables &,
               const eckit::Configuration &,
               const Parameters_ &,
               const oops::FieldSet3D &,
               const oops::FieldSet3D &);
  virtual ~ShadowLevels() = default;

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
