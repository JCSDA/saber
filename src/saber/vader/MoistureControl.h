/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class MoistureControlCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MoistureControlCovarianceParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::Parameter<int> mu_bins{"rht bins", "relative humidity bins", 30, this};
};

// -----------------------------------------------------------------------------

class MoistureControlParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistureControlParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<MoistureControlCovarianceParameters>
    moistureControlParams{"covariance data", this};
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({
        "qt",
        "mu",
        "potential_temperature",
        "virtual_potential_temperature"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"virtual_potential_temperature",
                          "mu"});
    const int modelLevels = outerVars.getLevels("qt");
    vars.addMetaData("virtual_potential_temperature", "levels", modelLevels);
    vars.addMetaData("mu", "levels", modelLevels);
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({"potential_temperature",
                          "qt"});
    for (const auto & var : vars.variables()) {
      vars.addMetaData(var, "levels", outerVars.getLevels(var));
    }
    return vars;
  }
};

// -----------------------------------------------------------------------------

class MoistureControl : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::MoistureControl";}

  typedef MoistureControlParameters Parameters_;

  MoistureControl(const oops::GeometryData &,
                  const oops::Variables &,
                  const eckit::Configuration &,
                  const Parameters_ &,
                  const oops::FieldSet3D &,
                  const oops::FieldSet3D &);
  virtual ~MoistureControl();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D & fset) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const atlas::FieldSet &,
                              const MoistureControlCovarianceParameters &);

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
