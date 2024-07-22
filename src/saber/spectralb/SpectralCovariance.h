/*
 * (C) Crown Copyright 2022-2024 Met Office
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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralCovarianceParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<spectralbCalibrationVertCovParameters>
    calibrationParams{"calibration", this};
  oops::OptionalParameter<spectralbReadParameters> readParams{"read", this};
  oops::Parameter<bool> skipVerticalConv{"skip vertical convolution",
    "Flag to skip vertical convolutions", false, this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class SpectralCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralCovariance";}

  typedef SpectralCovarianceParameters Parameters_;

  SpectralCovariance(const oops::GeometryData &,
                     const oops::Variables &,
                     const eckit::Configuration &,
                     const Parameters_ &,
                     const oops::FieldSet3D &,
                     const oops::FieldSet3D &);

  virtual ~SpectralCovariance() = default;

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void write() const override;

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;

  /// Active variables
  const oops::Variables activeVars_;

  /// netCDF configuration created for setting the netCDF header
  /// currently :: used in calibration mode
  /// TODO (Marek) extend to covariance reading to check yaml
  ///              covariance file consistency.
  eckit::LocalConfiguration netCDFConf_;

  /// Vertical Spectral Covariances
  atlas::FieldSet spectralVerticalCovariances_;

  /// Geometry data
  const oops::GeometryData & geometryData_;

  /// Spectral FunctionSpace
  const atlas::functionspace::Spectral specFunctionSpace_;

  /// Latest date time in sample data :: used in calibration mode
  const util::DateTime datetime_;

  /// Number of samples :: used in calibration
  std::size_t sample_size_;
};

}  // namespace spectralb
}  // namespace saber
