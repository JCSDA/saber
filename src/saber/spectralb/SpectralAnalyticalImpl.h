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

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralAnalyticalImplParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralAnalyticalImplParameters, SaberBlockParametersBase)

 public:
    /// Whether to preserve variance of processed increments
    oops::Parameter<bool> preservingVariance{"preserving variance", false, this};

    /// Define filter as the complement of the function
    oops::Parameter<bool> complementFilter{"complement filter", false, this};

    /// Filter specifications (Gaussian, boxcar function, triangle...)
    oops::Parameter<eckit::LocalConfiguration> function{"function",
                                                        eckit::LocalConfiguration(), this};

    oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralAnalyticalImpl {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralAnalyticalImpl";}

  typedef SpectralAnalyticalImplParameters Parameters_;

  SpectralAnalyticalImpl(const oops::GeometryData &,
                         const oops::Variables &,
                         const eckit::Configuration &,
                         const Parameters_ &,
                         const util::DateTime &,
                         const bool &);

  const oops::GeometryData & innerGeometryData() const {return innerGeometryData_;}
  const oops::Variables & innerVars() const {return innerVars_;}

  void randomize(oops::FieldSet3D &) const
    {throw eckit::Exception("Not implemented yet", Here());}
  void multiply(oops::FieldSet3D &) const;
  void multiplyAD(oops::FieldSet3D &) const;
  void leftInverseMultiply(oops::FieldSet3D &) const;

  size_t ctlVecSize() const
    {return ctlVecSize_;}

  // For inverse tests
  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                         const oops::Variables & innerVars) const;

  // For inverse tests
  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & outerVars) const;

  // Print
  void print(std::ostream &) const;

 private:
  /// Parameters
  Parameters_ params_;
  /// Valid time
  util::DateTime validTime_;
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
  /// Control vector size
  size_t ctlVecSize_;
};

}  // namespace spectralb
}  // namespace saber
