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

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/vader/MoistureControlParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
class MoistureControlSaberBlockParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MoistureControlSaberBlockParameters, SaberOuterBlockParametersBase)
 public:
  oops::RequiredParameter<moisturecontrolParameters>
    moisturecontrolParams{"covariance data", this};
};

// -----------------------------------------------------------------------------
// This saber block is here
//
// -----------------------------------------------------------------------------
class MoistureControlSaberBlock : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::MoistureControlSaberBlock";}

  typedef MoistureControlSaberBlockParameters Parameters_;

  MoistureControlSaberBlock(const eckit::mpi::Comm &,
         const oops::GeometryData &,
         const std::vector<size_t> &,
         const eckit::Configuration &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);
  virtual ~MoistureControlSaberBlock();

  const oops::GeometryData & inputGeometryData() const override {return inputGeometryData_;}
  const oops::Variables & inputVars() const override {return inputVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & inputGeometryData_;
  oops::Variables inputVars_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace saber
