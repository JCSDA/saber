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

#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/spectralb/SpectralAnalyticalImpl.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralAnalyticalCorrelation : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralAnalyticalCorrelation";}

  typedef SpectralAnalyticalImplParameters Parameters_;

  SpectralAnalyticalCorrelation(const oops::GeometryData &,
                                const oops::Variables &,
                                const eckit::Configuration &,
                                const Parameters_ &,
                                const oops::FieldSet3D &,
                                const oops::FieldSet3D &);

  virtual ~SpectralAnalyticalCorrelation() = default;

  void randomize(oops::FieldSet3D & fset) const override
    {impl_->randomize(fset);}
  void multiply(oops::FieldSet3D & fset) const override;

  size_t ctlVecSize() const override
    {return impl_->ctlVecSize();}
  void multiplySqrt(const atlas::Field &,
                    oops::FieldSet3D &,
                    const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &,
                      atlas::Field &,
                      const size_t &) const override;

 private:
  void print(std::ostream & os) const override
    {impl_->print(os);}

  /// Analytical implementation
  std::unique_ptr<spectralb::SpectralAnalyticalImpl> impl_;
};

}  // namespace spectralb
}  // namespace saber
