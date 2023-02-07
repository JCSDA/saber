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

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

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
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<MoistureControlCovarianceParameters>
    moistureControlParams{"covariance data", this};
};

// -----------------------------------------------------------------------------

class MoistureControl : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::MoistureControl";}

  typedef MoistureControlParameters Parameters_;

  MoistureControl(const oops::GeometryData &,
                  const std::vector<size_t> &,
                  const oops::Variables &,
                  const Parameters_ &,
                  const atlas::FieldSet &,
                  const atlas::FieldSet &,
                  const std::vector<atlas::FieldSet> &);
  virtual ~MoistureControl();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

atlas::FieldSet createMuStats(const atlas::FieldSet &,
                              const MoistureControlCovarianceParameters &);

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
