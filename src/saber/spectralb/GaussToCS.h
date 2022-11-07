/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/functionspace.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace interpolation {
class AtlasInterpWrapper;
}
}


namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------
class GaussToCSParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussToCSParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<oops::Variables> activeVariables{"active variables", this};
  // No parameters for now (in the future may add N as a parameter if it is possible
  // to use the one different from the one inferred from the gaussian grid
  oops::RequiredParameter<std::string> gaussGridUid{"gauss grid uid",
    "Gauss Grid UID", this};
};

// -----------------------------------------------------------------------------

class GaussToCS : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::GaussToCS";}

  typedef GaussToCSParameters Parameters_;

  GaussToCS(const oops::GeometryData &,
            const std::vector<size_t> &,
            const oops::Variables &,
            const Parameters_ &,
            const atlas::FieldSet &,
            const atlas::FieldSet &,
            const std::vector<atlas::FieldSet> &);
  virtual ~GaussToCS() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

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

  const oops::GeometryData innerGeometryData_;
};

}  // namespace spectralb;
}  // namespace saber
