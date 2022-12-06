/*
 * (C) Copyright 2022- UCAR
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/trans/Trans.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------
class SpectralToGaussParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralToGaussParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<oops::Variables> activeVariables{"active variables", this};
  // In the future may add N as a parameter if it is possible
  // to use the one different from the one inferred from the gaussian grid
};

// -----------------------------------------------------------------------------

class SpectralToGauss : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralToGauss";}

  typedef SpectralToGaussParameters Parameters_;

  SpectralToGauss(const oops::GeometryData &,
                  const std::vector<size_t> &,
                  const oops::Variables &,
                  const Parameters_ &,
                  const atlas::FieldSet &,
                  const atlas::FieldSet &,
                  const std::vector<atlas::FieldSet> &);
  virtual ~SpectralToGauss() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  void multiplyVectorFields(atlas::FieldSet &, atlas::FieldSet &) const;
  void multiplyVectorFieldsAD(atlas::FieldSet &, atlas::FieldSet &) const;
  void multiplyScalarFields(const atlas::FieldSet &, atlas::FieldSet &) const;
  void multiplyScalarFieldsAD(const atlas::FieldSet &, atlas::FieldSet &) const;

  Parameters_ params_;
  oops::Variables innerVars_;
  oops::Variables outerVars_;
  oops::Variables activeVars_;
  std::vector<std::size_t> activeVariableSizes_;

  /// Whether to convert to and from u/v.
  const bool useWindTransform;
  /// Gaussian (outer) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  /// Spectral (inner) functionspace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Trans object for gaussian-spectral transforms
  const atlas::trans::Trans trans_;

  const oops::GeometryData innerGeometryData_;
};

}  // namespace spectralb
}  // namespace saber
