/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/interpolation/Rescaling.h"

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------
class GaussToCSParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussToCSParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<oops::Variables> activeVariables{"active variables", this};
  // No parameters for now (in the future may add N as a parameter if it is possible
  // to use the one different from the one inferred from the gaussian grid
  oops::RequiredParameter<std::string> gaussGridUid{"gauss grid uid",
    "Gauss Grid UID", this};
  oops::Parameter<bool> initializeInverseInterpolation{"initialize inverse interpolator",
    true, this};
  oops::OptionalParameter<eckit::LocalConfiguration> interpolationRescaling{
    "rescaling", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------


struct CS2Gauss {
  // PointCloud FunctionSpace on Gauss grid with matching CubedSphere partitioner
  std::unique_ptr<const atlas::functionspace::PointCloud> matchingPtcldFspace;
  // PointCloud FunctionSpace on Gauss grid with Gauss grid partitioner
  std::unique_ptr<const atlas::functionspace::PointCloud> targetPtcldFspace;
  // Interpolation from CubedSpere NodeColumns to matchingPtcldFspace
  atlas::Interpolation interpolation;
  // Redistribution of MPI tasks from matchingPtcldFspace to targetPtcldFspace
  atlas::Redistribution redistribution;
};

// ------------------------------------------------------------------------------

class GaussToCS : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::GaussToCS";}

  typedef GaussToCSParameters Parameters_;

  GaussToCS(const oops::GeometryData &,
            const oops::Variables &,
            const eckit::Configuration &,
            const Parameters_ &,
            const oops::FieldSet3D &,
            const oops::FieldSet3D &);
  virtual ~GaussToCS() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                         const oops::Variables & innerVars) const override;

  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & outerVars) const override;

 private:
  void print(std::ostream &) const override;

  oops::Variables innerVars_;
  oops::Variables activeVars_;

  /// Cubed-sphere-dual NodeColumns space (outer) functionspace
  const atlas::functionspace::NodeColumns CSFunctionSpace_;

  /// Gaussian (grid)
  const atlas::StructuredGrid gaussGrid_;

  /// Gaussian (inner) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;

  /// Gaussian Partitioner
  const atlas::grid::Partitioner gaussPartitioner_;

  /// Cubed-sphere grid (destination grid)
  const atlas::Grid csgrid_;

  /// Interpolation Wrapper
  saber::interpolation::AtlasInterpWrapper interp_;

  /// Wrapper for inverse interpolation objects
  const CS2Gauss inverseInterpolation_;

  /// Optional rescaling weights
  const saber::interpolation::Rescaling rescaling_;

  const oops::GeometryData innerGeometryData_;
};

}  // namespace interpolation
}  // namespace saber
