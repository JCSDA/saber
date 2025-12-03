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

#include "saber/bifourier/BifourierTransformBase.h"
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
};

// -----------------------------------------------------------------------------

class BifourierBalanceCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierBalanceCalibrationParameters, oops::Parameters)

 public:
  // Use full recursive inverse formula to compute the regression
  oops::Parameter<bool> fullRecursiveInverse{"full recursive inverse", false, this};

  // Filtering scale (in total wavenumber unit)
  oops::Parameter<double> filteringScale{"filtering scale", 0.0, this};

  // Remaining variance fraction (between 0 and 1) in the auto-covariance inversion
  oops::Parameter<double> remainingVar{"remaining variance fraction", 1.0, this};

  // Old covariance input file
  oops::OptionalParameter<std::string> oldCovInputFile{"old covariance input file", this};

  // Half life
  oops::OptionalParameter<double> halfLife{"half life", this};

  // Cycle index
  oops::OptionalParameter<size_t> cycleIndex{"cycle index", this};

  // Sub-ensembles size
  oops::Parameter<size_t> subEnsSize{"sub-ensembles size", 0, this};
};

// -----------------------------------------------------------------------------

class BifourierBalanceWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierBalanceWriteParameters, oops::Parameters)

 public:
  // Output file
  oops::RequiredParameter<std::string> outputFile{"output file", this};

  // Write covariance flag
  oops::Parameter<bool> writeCovariance{"write covariance", false, this};
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

  // Calibration parameters
  oops::OptionalParameter<BifourierBalanceCalibrationParameters> calibration{"calibration", this};

  // Write parameters
  oops::OptionalParameter<BifourierBalanceWriteParameters> write{"write", this};

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

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  void write() const override;

 protected:
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Filtering length-scale
  const double Lf_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  const std::shared_ptr<BifourierTransformBase> trans_;

  // Ordered variables
  oops::Variables balVars_;

  // Number of active regression components
  size_t nCmp_;

  // Data
  atlas::FieldSet data_;

  // Interative counter
  size_t iterativeN_;

  // Private methods

  // Read covariance
  void readCovariance();

  // Compute regression
  void computeRegression(const std::vector<std::string> &,
                         const oops::Variable &);

  // Compute regressions from covariances
  void computeRegressionsFromCovariances();

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
