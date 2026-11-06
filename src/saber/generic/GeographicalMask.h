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

#include "oops/base/FieldSet4D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

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
  // Input mask ATLAS file
  oops::OptionalParameter<eckit::LocalConfiguration> inputMaskAtlasFileConf{
    "input mask atlas file", this};

  // Input mask model file
  oops::OptionalParameter<eckit::LocalConfiguration> inputMaskModelFileConf{
    "input mask model file", this};

  // Input weight ATLAS file
  oops::OptionalParameter<eckit::LocalConfiguration> inputWgtAtlasFileConf{
    "input weight atlas file", this};

  // Input weight model file
  oops::OptionalParameter<eckit::LocalConfiguration> inputWgtModelFileConf{
    "input weight model file", this};

  // Output weight ATLAS file
  oops::OptionalParameter<eckit::LocalConfiguration> outputWgtAtlasFileConf{
    "output weight atlas file", this};

  // Output weight model file
  oops::OptionalParameter<eckit::LocalConfiguration> outputWgtModelFileConf{
    "output weight model file", this};

  // Mask variable
  oops::OptionalParameter<std::string> maskVariable{"mask variable", this};

  // Number of levels of the mask
  oops::OptionalParameter<size_t> maskLevels{"mask levels", this};

  // Masks parameters vector
  oops::RequiredParameter<std::vector<MaskParameters>> masksParams{"masks", this};

  // Masks filtering
  oops::OptionalParameter<std::vector<OuterBlockParametersWrapper>> filterParams{
    "filter", this};
};

// -----------------------------------------------------------------------------

class GeographicalMaskParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GeographicalMaskParameters, BlockParametersBase)

 public:
  oops::OptionalParameter<GeographicalMaskParametersBase> read{"read", this};
  oops::OptionalParameter<GeographicalMaskParametersBase> calibration{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class GeographicalMask : public OuterBlockBase {
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

  void write() const override;

 private:
  const oops::GeometryData & innerGeometryData_;
  oops::Variables activeOuterVars_;
  ParametersBase_ params_;
  oops::Variables innerVars_;

  // Background and first guess for the filter constructor
  std::unique_ptr<oops::FieldSet4D> fset4dXb_;
  std::unique_ptr<oops::FieldSet4D> fset4dFg_;

  // Mask
  std::unique_ptr<oops::FieldSet3D> maskFset_;

  // Weight
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
