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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/bifourier/BifourierTransformBase.h"
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

  // Read covariance flag
  oops::Parameter<bool> readCovariance{"read covariance", false, this};
};

// -----------------------------------------------------------------------------

class ProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProfileParameters, oops::Parameters)

 public:
  // Variable name
  oops::RequiredParameter<std::string> variable{"variable", this};

  // Horizontal length-scales
  oops::RequiredParameter<std::vector<double>> Lh{"horizontal length-scales", this};

  // Vertical length-scale
  oops::RequiredParameter<double> Lv{"vertical length-scale", this};

  // Standard-deviation
  oops::OptionalParameter<std::vector<double>> stdDev{"standard-deviation", this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceCalibrationParameters, oops::Parameters)

 public:
  // Filtering scale (in total wavenumber unit)
  oops::Parameter<double> filteringScale{"filtering scale", 0.0, this};

  // Old covariance input file
  oops::OptionalParameter<std::string> oldCovInputFile{"old covariance input file", this};

  // Half life
  oops::OptionalParameter<double> halfLife{"half life", this};

  // Cycle index
  oops::OptionalParameter<size_t> cycleIndex{"cycle index", this};

  // Sub-ensembles size
  oops::Parameter<size_t> subEnsSize{"sub-ensembles size", 0, this};

  // User-defined vertical profile for each variable
  oops::Parameter<std::vector<ProfileParameters>> profiles{"profiles", {}, this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceWriteParameters, oops::Parameters)

 public:
  // Output file
  oops::RequiredParameter<std::string> outputFile{"output file", this};

  // Write covariance flag
  oops::Parameter<bool> writeCovariance{"write covariance", false, this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierCovarianceReadParameters> read{"read", this};

  // Calibration parameters
  oops::OptionalParameter<BifourierCovarianceCalibrationParameters> calibration{"calibration",
    this};

  // Standard-deviation inflation factor
  oops::Parameter<double> inflation{"inflation", 1.0, this};

  // Only correlation
  oops::Parameter<bool> correlation{"correlation", false, this};

  // Write parameters
  oops::OptionalParameter<BifourierCovarianceWriteParameters> write{"write", this};

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

  size_t ctlVecSize() const override
    {return trans_->ctlVecSize();}
  void randomCtlVec(atlas::Field & cv,
                    const size_t & offset) const override
    {trans_->randomCtlVec(cv, centralVars(), offset);}
  void multiplySqrt(const atlas::Field &,
                    oops::FieldSet3D &,
                    const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &,
                      atlas::Field &,
                      const size_t &) const override;

  void read() override;

  void directCalibration(const oops::FieldSets &) override;

  void iterativeCalibrationInit() override;
  void iterativeCalibrationUpdate(const oops::FieldSet3D &) override;
  void iterativeCalibrationFinal() override;

  void write() const override;

 protected:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Parameters
  Parameters_ params_;

  // Filtering length-scale
  const double Lf_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  const std::shared_ptr<BifourierTransformBase> trans_;

  // Data
  atlas::FieldSet data_;

  // Interative counter
  size_t iterativeN_;

  // Private methods

  // Read covariance
  void readCovariance();

  // Compute square-root
  void computeSquareRoot();

  // Compute covariance from correlation square-root and standard-deviation
  void computeCovariance(atlas::FieldSet &) const;

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
