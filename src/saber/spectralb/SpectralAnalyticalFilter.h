/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralAnalyticalFilterParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralAnalyticalFilterParameters, SaberBlockParametersBase)

 public:
    /// Whether to normalize as a localization function
    oops::Parameter<bool> normalizeFilterVariance{"normalize filter variance", false, this};

    /// Whether to preserve variance of processed increments
    /// Should only be used when "normalize filter variance" == false.
    oops::Parameter<bool> preservingVariance{"preserving variance", false, this};

    /// Define filter as the complement of the function
    oops::Parameter<bool> complementFilter{"complement filter", false, this};

    /// Filter specifications (Gaussian, boxcar function, triangle...)
    oops::Parameter<eckit::LocalConfiguration> function{"function",
                                                        eckit::LocalConfiguration(), this};

    oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralAnalyticalFilter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralAnalyticalFilter";}

  typedef SpectralAnalyticalFilterParameters Parameters_;

  SpectralAnalyticalFilter(const oops::GeometryData &,
                           const oops::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~SpectralAnalyticalFilter() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  // For inverse tests
  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                         const oops::Variables & innerVars) const override;

  // For inverse tests
  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & outerVars) const override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::Variables activeVars_;
  /// inner Geometry Data for next block
  const oops::GeometryData & innerGeometryData_;
  /// inner variables for next block
  const oops::Variables innerVars_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Filter in spectral space
  const std::vector<double> spectralFilter_;
};

}  // namespace spectralb
}  // namespace saber
