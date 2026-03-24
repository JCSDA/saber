/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/spectralb/SpectralAnalyticalImpl.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralAnalyticalFilter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralAnalyticalFilter";}

  typedef SpectralAnalyticalImplParameters Parameters_;

  SpectralAnalyticalFilter(const oops::GeometryData &,
                           const oops::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~SpectralAnalyticalFilter() = default;

  const oops::GeometryData & innerGeometryData() const override {return impl_->innerGeometryData();}
  const oops::Variables & innerVars() const override {return impl_->innerVars();}

  void multiply(oops::FieldSet3D & fset) const override
    {impl_->multiply(fset);}
  void multiplyAD(oops::FieldSet3D & fset) const override
    {impl_->multiplyAD(fset);}
  void leftInverseMultiply(oops::FieldSet3D & fset) const override
    {impl_->leftInverseMultiply(fset);}

  // For inverse tests
  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                         const oops::Variables & innerVars) const override
    {return impl_->generateInnerFieldSet(innerGeometryData, innerVars);}

  // For inverse tests
  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & outerVars) const override
    {return impl_->generateOuterFieldSet(outerGeometryData, outerVars);}

 private:
  void print(std::ostream & os) const override
    {impl_->print(os);}

  /// Analytical implementation
  std::unique_ptr<spectralb::SpectralAnalyticalImpl> impl_;
};

}  // namespace spectralb
}  // namespace saber
