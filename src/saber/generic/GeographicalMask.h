/*
 * (C) Copyright 2026 Meteorologisk Institutt
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
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

enum class RelationalOperator {
  LESS_THAN,
  LESS_THAN_OR_EQUAL,
  EQUAL,
  NOT_EQUAL,
  GREATER_THAN_OR_EQUAL,
  GREATER_THAN
};

struct RelationalOperatorParameterTraitsHelper {
  typedef RelationalOperator EnumType;
  static constexpr char enumTypeName[] = "RelationalOperator";
  static constexpr util::NamedEnumerator<RelationalOperator> namedValues[] = {
    { RelationalOperator::LESS_THAN, "<" },
    { RelationalOperator::LESS_THAN_OR_EQUAL, "<=" },
    { RelationalOperator::EQUAL, "==" },
    { RelationalOperator::NOT_EQUAL, "!=" },
    { RelationalOperator::GREATER_THAN_OR_EQUAL, ">=" },
    { RelationalOperator::GREATER_THAN, ">" }
  };
};

// -----------------------------------------------------------------------------

class MaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaskParameters, oops::Parameters)

 public:
  // Suffix
  oops::RequiredParameter<std::string> suffix{"suffix", this};

  // Relational operator
  oops::OptionalParameter<RelationalOperator> relationalOperator{"relational operator", this};

  // Threshold
  oops::OptionalParameter<double> threshold{"threshold", this};
};

// -----------------------------------------------------------------------------

class GeographicalMaskParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeographicalMaskParametersBase, oops::Parameters)

 public:
  // Input mask file
  oops::OptionalParameter<eckit::LocalConfiguration> inputMaskFileConf{
    "input mask file", this};

  // Input weight file
  oops::OptionalParameter<eckit::LocalConfiguration> inputWgtFileConf{
    "input weight file", this};

  // Output weight file
  oops::OptionalParameter<eckit::LocalConfiguration> outputWgtFileConf{
    "output weight file", this};

  // Mask variable
  oops::OptionalParameter<std::string> maskVariable{"mask variable", this};

  // Masks parameters vector
  oops::RequiredParameter<std::vector<MaskParameters>> masksParams{"masks", this};
};

// -----------------------------------------------------------------------------

class GeographicalMaskParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GeographicalMaskParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<GeographicalMaskParametersBase> read{"read", this};
  oops::OptionalParameter<GeographicalMaskParametersBase> calibration{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class GeographicalMask : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::generic::GeographicalMask";}

  typedef GeographicalMaskParameters     Parameters_;
  typedef GeographicalMaskParametersBase ParametersBase_;

  GeographicalMask(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~GeographicalMask() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return innerGeometryData_;}
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
  const oops::GeometryData & innerGeometryData_;
  oops::Variables activeOuterVars_;
  ParametersBase_ params_;
  oops::Variables innerVars_;

  // Mask
  eckit::LocalConfiguration readMaskConf_;
  std::unique_ptr<oops::FieldSet3D> maskFset_;

  // Weight
  eckit::LocalConfiguration readWgtConf_;
  eckit::LocalConfiguration writeWgtConf_;
  std::unique_ptr<oops::FieldSet3D> wgtFset_;

  void print(std::ostream &) const override;
};

}  // namespace generic
}  // namespace saber

// -----------------------------------------------------------------------------

namespace oops {
  template <>
  struct ParameterTraits<saber::generic::RelationalOperator> :
    public EnumParameterTraits<saber::generic::RelationalOperatorParameterTraitsHelper>{};
}  // namespace oops
