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
#include "oops/util/parameters/Parameters.h"

#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierBalanceReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierBalanceReadParameters, oops::Parameters)

 public:
  // Input file
  oops::RequiredParameter<std::string> inputFile{"input file", this};

  // Input file format ("netcdf", "arome legacy binary" or "arome legacy netcdf")
  oops::Parameter<std::string> inputFileFormat{"input file format", "netcdf", this};
};

// -----------------------------------------------------------------------------

class BifourierBalanceRowParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierBalanceRowParameters, oops::Parameters)

 public:
  // Output variable
  oops::RequiredParameter<std::string> outputVar{"output variable", this};

  // Input variables
  oops::Parameter<std::vector<std::string>> inputVars{"input variables", {}, this};
};

// -----------------------------------------------------------------------------

class BifourierBalanceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierBalanceParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierBalanceReadParameters> read{"read", this};

  // Output file
  oops::OptionalParameter<std::string> outputFile{"output file", this};

  // Rows
  oops::RequiredParameter<std::vector<BifourierBalanceRowParameters>>
    rows{"rows", this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierBalance : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierBalance";}

  typedef BifourierBalanceParameters Parameters_;

  BifourierBalance(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~BifourierBalance() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return innerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override;

  void write() const override;

 private:
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  const oops::Variables & innerVars_;

  // Parameters
  BifourierBalanceParameters params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransform> trans_;

  // Ordered variables
  oops::Variables balVars_;

  // Number of active regression components
  size_t nCmp_;

  // Data
  atlas::FieldSet data_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
