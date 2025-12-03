/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>

#include "saber/bifourier/BifourierCovariance.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierAromeCovarianceReadParameters : public BifourierCovarianceReadParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeCovarianceReadParameters,
    BifourierCovarianceReadParameters)

 public:
  // Input file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> inputFileFormat{"input file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BifourierAromeCovarianceWriteParameters : public BifourierCovarianceWriteParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeCovarianceWriteParameters,
    BifourierCovarianceWriteParameters)

 public:
  // Output file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> outputFileFormat{"output file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BifourierAromeCovarianceParameters : public BifourierCovarianceParameters {
  OOPS_CONCRETE_PARAMETERS(BifourierAromeCovarianceParameters, BifourierCovarianceParameters)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierAromeCovarianceReadParameters> read{"read", this};

  // Write parameters
  oops::OptionalParameter<BifourierAromeCovarianceWriteParameters> write{"write", this};

  // REDNMC factor
  oops::Parameter<double> rednmc{"rednmc", std::sqrt(0.5), this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};


// -----------------------------------------------------------------------------

class BifourierAromeCovariance : public BifourierCovariance {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierAromeCovariance";}

  typedef BifourierAromeCovarianceParameters Parameters_;

  BifourierAromeCovariance(const oops::GeometryData & gdata,
                           const oops::Variables & activeVars,
                           const eckit::Configuration & covarConf,
                           const Parameters_ & params,
                           const oops::FieldSet3D & xb,
                           const oops::FieldSet3D & fg) :
    BifourierCovariance(gdata, activeVars, covarConf, params, xb, fg), params_(params) {}

  void read();

  void write() const;

 private:
  // Parameters
  Parameters_ params_;

  // Private methods
  double aromeWeight(const size_t &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
