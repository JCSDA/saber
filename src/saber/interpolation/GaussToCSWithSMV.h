/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>

#include "atlas/interpolation/Interpolation.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/interpolation/Rescaling.h"
#include "saber/interpolation/SMVInterpWrapper.h"

namespace saber {
namespace interpolation {

// WARNING: This Saber block currently only works if run on "Rubik's Cube" MPI
// ranks, e.g. 1, 6, 9x6, 25x6
// This is due to an issue in atlas when constructing a Gauss grid with a
// CubedSphere distribution - this only works if the region around each pole is
// contained within a single partition.


// -----------------------------------------------------------------------------
class GaussToCSWithSMVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussToCSWithSMVParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<oops::Variables> activeVariables{"active variables",
                                                           this};
  // No parameters for now (in the future may add N as a parameter if it is
  // possible to use the one different from the one inferred from the gaussian
  // grid
  oops::RequiredParameter<std::string> gaussGridUid{"gauss grid uid",
                                                    "Gauss Grid UID", this};
  // Interpolation type, defaults to "structured-bilinear" or
  // "unstructured-bilinear-lonlat"
  oops::Parameter<std::string> interpType{"interpolation type", "", this};
  oops::Parameter<bool> initializeInverseInterpolation{
      "initialize inverse interpolator", true, this};
  oops::Parameter<bool> normalizeInverseInterpolation{
      "normalize inverse interpolator", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> interpolationRescaling{
      "rescaling", this};
  // switch to treat winds as vectors, defaults to true
  oops::Parameter<bool> interpolateWindAsVector{
      "interpolate wind as vector", true, this};
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables();
  }
};

// ------------------------------------------------------------------------------

class GaussToCSWithSMV : public SaberOuterBlockBase {
 public:
  static const std::string classname() { return "saber::spectralb::GaussToCSWithSMV"; }

  typedef GaussToCSWithSMVParameters Parameters_;

  GaussToCSWithSMV(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration&,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~GaussToCSWithSMV() = default;

  const oops::GeometryData &innerGeometryData() const override {
    return innerGeometryData_;
  }
  const oops::Variables &innerVars() const override { return innerVars_; }

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &fset) const override {
    inverseMultiply(fset);
  }
  void rightInverseMultiply(oops::FieldSet3D &fset) const override {
    inverseMultiply(fset);
  }
  void directCalibration(const oops::FieldSets &) override;

  oops::FieldSet3D generateInnerFieldSet(
      const oops::GeometryData &innerGeometryData,
      const oops::Variables &innerVars) const override;

  oops::FieldSet3D generateOuterFieldSet(
      const oops::GeometryData &outerGeometryData,
      const oops::Variables &outerVars) const override;

 private:
  void print(std::ostream &) const override;
  void inverseMultiply(oops::FieldSet3D &) const;

  oops::Variables innerVars_;
  oops::Variables activeVars_;

  /// Cubed-sphere-dual NodeColumns space (outer) functionspace
  const atlas::functionspace::NodeColumns CSFunctionSpace_;

  /// Gaussian (grid)
  const atlas::StructuredGrid gaussGrid_;

  /// Interpolation type
  std::string interpType_;

  /// Gaussian Partitioner
  const atlas::grid::Partitioner gaussPartitioner_;

  /// Gaussian (inner) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;

  /// Cubed-sphere grid (destination grid)
  const atlas::Grid csgrid_;

  /// Switch to vector wind interpolation
  const bool includingVectorInterpolation_;

  /// Interpolation Wrapper
  saber::interpolation::AtlasInterpWrapper interp_;

  /// Wrapper for inverse interpolation objects
  std::optional<SMVInterpWrapper> inverseInterpolation_;

  /// Optional rescaling weights
  const saber::interpolation::Rescaling rescaling_;

  const oops::GeometryData innerGeometryData_;
};

}  // namespace interpolation
}  // namespace saber
