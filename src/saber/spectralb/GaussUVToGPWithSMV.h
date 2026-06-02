/*
 * (C) Crown Copyright 2026- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

// WARNING: This Saber block currently only works if run on "Rubik's Cube" MPI
// ranks, e.g. 1, 6, 9x6, 25x6
// This is due to an issue in atlas when constructing a Gauss grid with a
// CubedSphere distribution - this only works if the region around each pole is
// contained within a single partition.

class GaussUVToGPWithSMVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussUVToGPWithSMVParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<std::string> modelGridName{"model grid name", this};
  oops::OptionalParameter<std::string> gaussState{"gauss state", this};
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({std::vector<std::string>{
        "eastward_wind", "geostrophic_pressure_levels_minus_one",
        "northward_wind"}});
  }

  oops::Variables activeInnerVars(
      const oops::Variables &outerVars) const override {
    const int modelLevels = outerVars["eastward_wind"].getLevels();
    eckit::LocalConfiguration conf;
    conf.set("levels", modelLevels);
    oops::Variables vars;
    vars.push_back(oops::Variable{"eastward_wind", conf});
    vars.push_back(oops::Variable{"northward_wind", conf});
    return vars;
  }

  oops::Variables activeOuterVars(
      const oops::Variables &outerVars) const override {
    oops::Variables vars{{outerVars["geostrophic_pressure_levels_minus_one"]}};
    return vars;
  }

  oops::Parameter<bool> normalizeInverseInterpolation{
      "normalize inverse interpolator", false, this};
};

// -----------------------------------------------------------------------------
/// \brief saber block that converts zonal and meridional wind to geostrophic
///        pressure on a Gaussian latitude mesh.

class GaussUVToGPWithSMV : public SaberOuterBlockBase {
 public:
  static const std::string classname() {
    return "saber::spectralb::GaussUVToGP";
  }

  typedef GaussUVToGPWithSMVParameters Parameters_;

  GaussUVToGPWithSMV(const oops::GeometryData &, const oops::Variables &,
              const eckit::Configuration &, const Parameters_ &,
              const oops::FieldSet3D &, const oops::FieldSet3D &);

  virtual ~GaussUVToGPWithSMV() = default;

  const oops::GeometryData &innerGeometryData() const override {
    return innerGeometryData_;
  }
  const oops::Variables &innerVars() const override { return innerVars_; }

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  Parameters_ params_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;

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
