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

class MoistureControlReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MoistureControlReadParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::Parameter<int> mu_bins{"rht bins", "relative humidity bins", 30, this};
};

// -----------------------------------------------------------------------------

class MoistureControlParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistureControlParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<MoistureControlReadParameters> readParams{"read", this};

  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({std::vector<std::string>{
        "qt",
        "mu",
        "air_potential_temperature",
        "virtual_potential_temperature"}});
  }

  const oops::Variables mandatoryStateVars() const override {
    return oops::Variables({"qt", "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
                            "air_potential_temperature", "dimensionless_exner_function",
                            "dlsvpdT", "qsat", "rht"});
  }

  oops::Variables activeInnerVars(const oops::Variables& outerVars) const override {
    const int modelLevels = outerVars["qt"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back({"virtual_potential_temperature", conf});
    vars.push_back({"mu", conf});
    return vars;
  }

  oops::Variables activeOuterVars(const oops::Variables& outerVars) const override {
    oops::Variables vars({outerVars["air_potential_temperature"],
                          outerVars["qt"]});
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
  void read() override;
  void directCalibration(const oops::FieldSets & fset) override;
  void write() const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;
  const std::size_t nlevs_;
  Parameters_ params_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const size_t &,
                              const atlas::FieldSet &,
                              const MoistureControlReadParameters &);

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
