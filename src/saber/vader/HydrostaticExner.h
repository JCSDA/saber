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

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class HydrostaticExnerCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerCovarianceParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::RequiredParameter<int> covariance_nlat{"number of covariance latitude rings", this};
  oops::Parameter<int> gp_regression_bins{"gp regression bins", "gP regression bins", 18, this};
};

// -----------------------------------------------------------------------------

class HydrostaticExnerParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerParameters, SaberOuterBlockParametersBase)
 public:
  oops::RequiredParameter<HydrostaticExnerCovarianceParameters>
    hydrostaticexnerParams{"covariance data", this};
};

// -----------------------------------------------------------------------------
// This saber block is here to do 3 jobs
// 1) the vertical regression on geostrophic pressure
// 2) summing the result with unbalanced pressure to create hydrostatic_pressure
// 3) converting hydrostatic pressure to exner pressure.
// -----------------------------------------------------------------------------

class HydrostaticExner : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydrostaticExner";}

  typedef HydrostaticExnerParameters Parameters_;

  HydrostaticExner(const oops::GeometryData &,
                   const std::vector<size_t> &,
                   const oops::Variables &,
                   const Parameters_ &,
                   const atlas::FieldSet &,
                   const atlas::FieldSet &,
                   const std::vector<atlas::FieldSet> &);
  virtual ~HydrostaticExner();

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

atlas::FieldSet createGpRegressionStats(const atlas::FunctionSpace &,
                                        const atlas::FieldSet &,
                                        const std::vector<size_t> &,
                                        const oops::Variables &,
                                        const HydrostaticExnerCovarianceParameters &);

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
