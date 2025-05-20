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

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSets.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variable.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierCovarianceReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceReadParameters, oops::Parameters)

 public:
  // Input file
  oops::RequiredParameter<std::string> inputFile{"input file", this};

  // Input file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> inputFileFormat{"input file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierCovarianceReadParameters> read{"read", this};

  // Standard-deviation inflation factor
  oops::Parameter<double> inflation{"inflation", 1.0, this};

  // Only correlation
  oops::Parameter<bool> correlation{"correlation", false, this};

  // Output file
  oops::OptionalParameter<std::string> outputFile{"output file", this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierCovariance";}

  typedef BifourierCovarianceParameters Parameters_;

  BifourierCovariance(const oops::GeometryData &,
                      const oops::Variables &,
                      const eckit::Configuration &,
                      const BifourierCovarianceParameters &,
                      const oops::FieldSet3D &,
                      const oops::FieldSet3D &);
  virtual ~BifourierCovariance();

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  size_t ctlVecSize() const override
    {return trans_->ctlVecSize();}

  void multiplySqrt(const atlas::Field &,
                    oops::FieldSet3D &,
                    const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &,
                      atlas::Field &,
                      const size_t &) const override;

  void read() override;

  void write() const override;

 private:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Active variables
  const oops::Variables activeVars_;

  // Parameters
  Parameters_ params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransform> trans_;

  // Data
  atlas::FieldSet data_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
