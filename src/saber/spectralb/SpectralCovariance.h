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

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbParameters> readParams{"read", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralCovariance";}

  typedef SpectralCovarianceParameters Parameters_;

  SpectralCovariance(const oops::GeometryData &,
                     const std::vector<size_t> &,
                     const oops::Variables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const atlas::FieldSet &,
                     const atlas::FieldSet &,
                     const size_t &);

  virtual ~SpectralCovariance() = default;

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  void read() override;

 private:
  void testUUtConsistency(const double &,
                          const double & adjointTolerance = 1.0e-12) const;
  void multiplyUMatrix(atlas::FieldSet &) const;
  void multiplyUMatrixAD(atlas::FieldSet &) const;
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Variables sizes
  const std::vector<size_t> variableSizes_;
  /// Active variables
  const oops::Variables activeVars_;
  /// Option to use vertical covariances or correlations
  bool variance_opt_;
  /// Covariance statistics
  // Note: only need vertical covariances or correlations from this;
  // probably can be gotten in the ctor and saved here instead of cs_
  std::unique_ptr<CovStat_ErrorCov> cs_;
  /// Geometry data
  const oops::GeometryData & geometryData_;
  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Time rank
  size_t timeRank_;
};

}  // namespace spectralb
}  // namespace saber
