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

#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class GroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupParameters, oops::Parameters)

 public:
  // Group suffix
  oops::Parameter<std::string> suffix{"suffix", "_shadowLevels", this};

  // Group name
  oops::OptionalParameter<std::string> name{"group name", this};

  // Group variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};

  // Number of shadow levels
  oops::RequiredParameter<size_t> nz{"number of shadow levels", this};

  // Vertical length-scale
  oops::OptionalParameter<double> rv{"vertical length-scale", this};

  // Component weight
  oops::Parameter<double> cmpWgt{"component weight", 1.0, this};
};

// -----------------------------------------------------------------------------

class ShadowLevelsParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ShadowLevelsParametersBase, oops::Parameters)

 public:
  // Groups of variables
  oops::RequiredParameter<std::vector<GroupParameters>> groups{"groups", this};

  // Input weight file
  oops::OptionalParameter<eckit::LocalConfiguration> inputWgtFileConf{
    "input weight file", this};

  // Output weight file
  oops::OptionalParameter<eckit::LocalConfiguration> outputWgtFileConf{
    "output weight file", this};
};

// -----------------------------------------------------------------------------

class ShadowLevelsParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ShadowLevelsParameters, BlockParametersBase)

 public:
  oops::OptionalParameter<ShadowLevelsParametersBase> read{"read", this};
  oops::OptionalParameter<ShadowLevelsParametersBase> calibration{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class ShadowLevels : public OuterBlockBase {
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
  oops::Variables activeOuterVars_;
  ParametersBase_ params_;
  eckit::LocalConfiguration fieldsMetaData_;
  oops::Variables innerVars_;

  // Groups of variables
  struct Group {
    std::string suffix_;
    std::string name_;
    std::string varInModelFile_;
    std::vector<std::string> variables_;
    size_t nz_;
    double rv_;
    double cmpWgtSqrt_;
  };
  std::vector<Group> groups_;

  // Factor for the GC99 function
  const double autoConvolFactor_ = 0.52;
  const double rv2L_ = 1.0/3.53;

  // Weight
  eckit::LocalConfiguration readConf_;
  eckit::LocalConfiguration writeConf_;
  std::unique_ptr<oops::FieldSet3D> wgtFset_;

  void print(std::ostream &) const override;
};

}  // namespace generic
}  // namespace saber
