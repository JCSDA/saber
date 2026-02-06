/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>
#include <vector>

#include "saber/bifourier/BifourierBalance.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierAromeBalanceReadParameters : public BifourierBalanceReadParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeBalanceReadParameters, BifourierBalanceReadParameters)

 public:
  // Input file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> inputFileFormat{"input file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BifourierAromeBalanceWriteParameters : public BifourierBalanceWriteParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeBalanceWriteParameters, BifourierBalanceWriteParameters)

 public:
  // Output file
  oops::RequiredParameter<std::string> outputFile{"output file", this};

  // Output file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> outputFileFormat{"output file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BalancedAirPressureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BalancedAirPressureParameters, oops::Parameters)

 public:
  // Zonal wavenumbers size
  oops::RequiredParameter<size_t> M{"zonal truncation", this};

  // Meridional wavenumbers size
  oops::RequiredParameter<size_t> N{"meridional truncation", this};

  // Mean latitude
  oops::RequiredParameter<double> meanLat{"mean latitude", this};
};

// -----------------------------------------------------------------------------

class BifourierAromeBalanceParameters : public BifourierBalanceParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeBalanceParameters, BifourierBalanceParameters)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierAromeBalanceReadParameters> read{"read", this};

  // Write parameters
  oops::OptionalParameter<BifourierAromeBalanceWriteParameters> write{"write", this};

  // Explicit balanced air pressure parameters
  oops::OptionalParameter<BalancedAirPressureParameters>
    explicitPb{"explicit balanced air pressure parameters", this};

  // Balanced air pressure parameters from grid
  oops::Parameter<bool> pbFromTrans{"balanced air pressure parameters from grid", false, this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables(
    std::vector<std::string>({
    "air_upward_absolute_vorticity",
    "air_temperature",
    "log_of_air_pressure_at_surface",
    "air_temperature_and_log_of_air_pressure_at_surface"}));}
};

// -----------------------------------------------------------------------------

class BifourierAromeBalance : public BifourierBalance {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierAromeBalance";}

  typedef BifourierAromeBalanceParameters Parameters_;

  BifourierAromeBalance(const oops::GeometryData &,
                        const oops::Variables &,
                        const eckit::Configuration &,
                        const Parameters_ &,
                        const oops::FieldSet3D &,
                        const oops::FieldSet3D &);

  const oops::Variables & innerVars() const override
    {return aromeInnerVars_;}

  void multiply(oops::FieldSet3D &) const;
  void multiplyAD(oops::FieldSet3D &) const;
  void leftInverseMultiply(oops::FieldSet3D &) const;

  void read();

  void directCalibration(const oops::FieldSets &);

  void iterativeCalibrationUpdate(const oops::FieldSet3D &);

  void write() const;

 private:
  // Parameters
  BifourierAromeBalanceParameters params_;

  // Number of levels
  size_t nz_;

  // Vorticity to balanced pressure factor
  std::vector<double> fact1_;

  // AROME balance inner variables
  oops::Variables aromeInnerVars_;

  // Private methods

  // Generic inner variables
  oops::Variables genericInnerVars(const oops::Variables &);

  // Vorticity to balanced pressure
  void vorToPb(oops::FieldSet3D &) const;

  // Vorticity to balanced pressure, adjoint
  void vorToPbAD(oops::FieldSet3D &) const;

  // Vorticity to balanced pressure, left inverse
  void vorToPbLeftInverse(oops::FieldSet3D &) const;

  // Remove balanced pressure
  void removePb(oops::FieldSet3D &) const;

  // Remove balanced pressure, adjoint
  void removePbAD(oops::FieldSet3D &) const;

  // Remove balanced pressure, left inverse
  void removePbLeftInverse(oops::FieldSet3D &) const;

  // Split TPs
  void splitTPs(oops::FieldSet3D &) const;

  // Gather TPs
  void gatherTPs(oops::FieldSet3D &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
