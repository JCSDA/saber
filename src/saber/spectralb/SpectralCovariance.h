/*
 * (C) Crown Copyright 2022 Met Office
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"
#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralCovarianceParameters : public SaberCentralBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralCovarianceParameters, SaberCentralBlockParametersBase)

 public:
  oops::RequiredParameter<spectralbParameters> spectralbParams{"spectralb", this};
};

// -----------------------------------------------------------------------------

class SpectralCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralCovariance";}

  typedef SpectralCovarianceParameters Parameters_;

  SpectralCovariance(const oops::GeometryData &,
                     const std::vector<size_t> &,
                     const oops::Variables &,
                     const Parameters_ &,
                     const atlas::FieldSet &,
                     const atlas::FieldSet &,
                     const std::vector<atlas::FieldSet> &);

  virtual ~SpectralCovariance() = default;

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;

  /// Active variables
  oops::Variables activeVars_;
  /// Option to use vertical covariances or correlations
  bool variance_opt_;
  /// Covariance statistics
  // Note: only need vertical covariances or correlations from this;
  // probably can be gotten in the ctor and saved here instead of cs_
  std::unique_ptr<CovStat_ErrorCov> cs_;
  // Spectral FunctionSpace
  atlas::functionspace::Spectral specFunctionSpace_;
};

}  // namespace spectralb
}  // namespace saber
