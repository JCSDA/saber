/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {


class VariableGroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableGroupParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> groupVariableName{"group variable name", this};
  oops::RequiredParameter<oops::Variables> groupComponents{"group components", this};
};

// -----------------------------------------------------------------------------

class DuplicateVariablesParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DuplicateVariablesParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::vector<VariableGroupParameters>>
   variableGroupParameters{"variable groupings", this};
  /// Variables groups
  oops::Variables mandatoryActiveVars() const override {
    oops::Variables vars;
    std::vector<VariableGroupParameters> gp = variableGroupParameters.value();
    for (std::size_t i = 0; i < gp.size(); ++i) {
      std::string s = gp[i].groupVariableName.value();
      oops::Variables v = gp[i].groupComponents.value();
      vars.push_back(s);
      vars += v;
    }
    return vars;}
};


// -----------------------------------------------------------------------------

class DuplicateVariables : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::DuplicateVariables";}

  typedef DuplicateVariablesParameters Parameters_;

  DuplicateVariables(const oops::GeometryData &,
                     const oops::Variables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const oops::FieldSet3D &,
                     const oops::FieldSet3D &);

  virtual ~DuplicateVariables() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  std::vector<VariableGroupParameters> groups_;
  oops::Variables outerVars_;
  oops::Variables activeVars_;
  oops::Variables innerVars_;
  const oops::GeometryData & innerGeometryData_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
