/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "saber/bifourier/BifourierTransformBase.h"
#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/bifourier/BiperiodizationImpl.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierSpectralVorDivToGridWindParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierSpectralVorDivToGridWindParameters, SaberBlockParametersBase)

 public:
  // Backward mode
  oops::Parameter<bool> backward{"backward mode", false, this};

  // Dipole test
  oops::Parameter<bool> dipoleTest{"dipole test", false, this};

  // Biperiodization parameters
  oops::OptionalParameter<BiperiodizationImplParameters> biperParams{"biperiodization",
    this};

  // Outer spherical winds
  oops::Parameter<bool> outerSphericalWinds{"outer spherical winds", false, this};

  // Transform parameters
  oops::Parameter<BifourierTransformParameters> transform{"transform",
    BifourierTransformParameters(), this};

  oops::Variables mandatoryActiveVars() const override {
    oops::Variables mandatoryActiveVars;
    mandatoryActiveVars.push_back("air_upward_absolute_vorticity");
    mandatoryActiveVars.push_back("air_horizontal_divergence");
    if (outerSphericalWinds) {
      mandatoryActiveVars.push_back("eastward_wind");
      mandatoryActiveVars.push_back("northward_wind");
    } else {
      mandatoryActiveVars.push_back("geographical_x_wind");
      mandatoryActiveVars.push_back("geographical_y_wind");
    }
    return mandatoryActiveVars;
  }
};

// -----------------------------------------------------------------------------

class BifourierSpectralVorDivToGridWind : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierSpectralVorDivToGridWind";}

  typedef BifourierSpectralVorDivToGridWindParameters Parameters_;

  BifourierSpectralVorDivToGridWind(const oops::GeometryData &,
                                    const oops::Variables &,
                                    const eckit::Configuration &,
                                    const Parameters_ &,
                                    const oops::FieldSet3D &,
                                    const oops::FieldSet3D &);
  virtual ~BifourierSpectralVorDivToGridWind() = default;

  const oops::GeometryData & innerGeometryData() const override;
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override
    {}

 private:
  // Inner geometry data
  std::unique_ptr<oops::GeometryData> innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransformBase> trans_;

  // FFT backend
  const std::string fftBackend_;

  // Number of levels
  size_t nz_;

  // Index js for (jk, jl) = (0, 0)
  int jsZero_;

  // Map factor and Jacobian coefficients FieldSet
  oops::FieldSet3D data_;

  // Private methods

  // Forward application
  void forward(oops::FieldSet3D &,
               const oops::Variables &) const;

  // Forward application, adjoint
  void forwardAD(oops::FieldSet3D &,
                 const oops::Variables &) const;

  // Backward application
  void backward(oops::FieldSet3D &,
                const oops::Variables &) const;

  // Backward application, adjoint
  void backwardAD(oops::FieldSet3D &,
                  const oops::Variables &) const;

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
