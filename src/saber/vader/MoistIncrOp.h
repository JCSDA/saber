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

#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {


// -----------------------------------------------------------------------------

class MoistIncrOpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistIncrOpParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> svp_file{"saturation vapour pressure file", this};
  oops::RequiredParameter<std::string> mio_file{"moisture incrementing operator file", this};
};

// -----------------------------------------------------------------------------

class MoistIncrOp : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::AirTemperature";}

  typedef MoistIncrOpParameters Parameters_;

  MoistIncrOp(const oops::GeometryData &,
              const std::vector<size_t> &,
              const oops::Variables &,
              const Parameters_ &,
              const atlas::FieldSet &,
              const atlas::FieldSet &,
              const std::vector<atlas::FieldSet> &);
  virtual ~MoistIncrOp();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
