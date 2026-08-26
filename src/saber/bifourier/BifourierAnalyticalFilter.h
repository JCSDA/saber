/*
 * (C) Crown Copyright 2023 Met Office
 * (C) Copyright 2025 Meteorologisk Institutt
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
#include "oops/util/parameters/OptionalParameter.h"

#include "saber/bifourier/BifourierTransformBase.h"
#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierAnalyticalFilterParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierAnalyticalFilterParameters, SaberBlockParametersBase)

 public:
    // Minimum waveband (optional)
    oops::OptionalParameter<double> wavebandMin{"waveband min", this};

    // Peak waveband
    oops::RequiredParameter<double> wavebandPeak{"waveband peak", this};

    // Maximum waveband (optional)
    oops::OptionalParameter<double> wavebandMax{"waveband max", this};

    // Inverse mode
    oops::Parameter<bool> inverseMode{"inverse mode", false, this};

    oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierAnalyticalFilter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::bifourier::BifourierAnalyticalFilter";}

  typedef BifourierAnalyticalFilterParameters Parameters_;

  BifourierAnalyticalFilter(const oops::GeometryData &,
                            const oops::Variables &,
                            const eckit::Configuration &,
                            const Parameters_ &,
                            const oops::FieldSet3D &,
                            const oops::FieldSet3D &);

  virtual ~BifourierAnalyticalFilter() = default;

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
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  oops::Variables innerVars_;

  /// Active variables
  const oops::Variables activeVars_;

  // Parameters
  Parameters_ params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  const std::shared_ptr<BifourierTransformBase> trans_;

  /// Filter in spectral space
  const std::vector<double> spectralFilter_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

}  // namespace bifourier
}  // namespace saber
