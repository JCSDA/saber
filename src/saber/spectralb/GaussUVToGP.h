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
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class GaussUVToGPParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPParameters, SaberBlockParametersBase)
 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const {return oops::Variables({
    "eastward_wind",
    "geostrophic_pressure_levels_minus_one",
    "northward_wind"});}
};

// -----------------------------------------------------------------------------


class GaussUVToGP : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::GaussUVToGP";}

  typedef GaussUVToGPParameters Parameters_;

  GaussUVToGP(const oops::GeometryData &,
              const std::vector<std::size_t> &,
              const oops::Variables &,
              const Parameters_ &,
              const atlas::FieldSet &,
              const atlas::FieldSet &,
              const std::vector<atlas::FieldSet> &);

  virtual ~GaussUVToGP() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;

  Parameters_ params_;
  oops::Variables innerVars_;
  oops::Variables outerVars_;
  std::vector<std::size_t> activeVariableSizes_;
  std::string modelGridName_;
  std::string gaussStateName_;

  /// Gaussian (outer) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  /// Spectral (inner) functionspace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Trans object for gaussian-spectral transforms
  const atlas::trans::Trans trans_;
  const oops::GeometryData innerGeometryData_;
  const atlas::FieldSet augmentedState_;
};

}  // namespace spectralb
}  // namespace saber
