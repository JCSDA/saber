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

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierCovarianceImplReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceImplReadParameters, oops::Parameters)

 public:
  // Input file
  oops::RequiredParameter<std::string> inputFile{"input file", this};

  // Input file from balance operator
  oops::Parameter<bool> inputFileFromBalance{"input file from balance", false, this};
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
  oops::Parameter<double> Lv{"vertical length-scale", 0.0, this};

  // Standard-deviation
  oops::OptionalParameter<std::vector<double>> stdDev{"standard-deviation", this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceImplCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceImplCalibrationParameters, oops::Parameters)

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
};

// -----------------------------------------------------------------------------

class BifourierCovarianceImplWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceImplWriteParameters, oops::Parameters)

 public:
  // Output file
  oops::RequiredParameter<std::string> outputFile{"output file", this};

  // Write covariance flag
  oops::Parameter<bool> writeCovariance{"write covariance", false, this};
};

// -----------------------------------------------------------------------------

class BifourierCovarianceImplParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierCovarianceImplParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<BifourierCovarianceImplReadParameters> read{"read", this};

  // Calibration parameters
  oops::OptionalParameter<BifourierCovarianceImplCalibrationParameters> calibration{"calibration",
    this};

  // User-defined vertical profile for each variable
  oops::OptionalParameter<std::vector<ProfileParameters>> profiles{"profiles", this};

  // Standard-deviation inflation factor
  oops::Parameter<double> inflation{"inflation", 1.0, this};

  // Only correlation
  oops::Parameter<bool> correlation{"correlation", false, this};

  // Write parameters
  oops::OptionalParameter<BifourierCovarianceImplWriteParameters> write{"write", this};

  oops::Variables mandatoryActiveVars() const
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierCovarianceImpl {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierCovarianceImpl";}

  typedef BifourierCovarianceImplParameters Parameters_;

  BifourierCovarianceImpl(const oops::GeometryData &,
                          const oops::Variables &,
                          const eckit::Configuration &,
                          const BifourierCovarianceImplParameters &,
                          const oops::FieldSet3D &,
                          const oops::FieldSet3D &);

  // Central and outer blocks methods
  const oops::GeometryData & innerGeometryData() const
    {return geometryData_;}
  const oops::Variables & innerVars() const
    {return vars_;}

  size_t ctlVecSize() const
    {return trans_->ctlVecSize();}
  void randomCtlVec(atlas::Field & cv,
                    const size_t & offset) const
    {trans_->randomCtlVec(cv, vars_, offset);}
  void multiplySqrt(oops::FieldSet3D &) const;
  void multiplySqrt(const atlas::Field &,
                    oops::FieldSet3D &,
                    const size_t &) const;
  void multiplySqrtAD(oops::FieldSet3D &) const;
  void multiplySqrtAD(const oops::FieldSet3D &,
                      atlas::Field &,
                      const size_t &) const;
  void leftInverseMultiply(oops::FieldSet3D &) const;

  void read();

  void directCalibration(const oops::FieldSets &);

  void iterativeCalibrationInit();
  void iterativeCalibrationUpdate(const oops::FieldSet3D &);
  void iterativeCalibrationFinal();

  void write() const;

  void print(std::ostream &) const;

  // Specific accessors
  const eckit::mpi::Comm & comm() const
    {return comm_;}
  atlas::FieldSet & data()
    {return data_;}
  const std::shared_ptr<BifourierTransformBase> & trans() const
    {return trans_;}

  // Compute square-root
  void computeSquareRoot();

  // Compute covariance from correlation square-root and standard-deviation
  void computeCovariance(atlas::FieldSet &) const;

 protected:
  // Geometry data
  const oops::GeometryData & geometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Variables
  const oops::Variables vars_;

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
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
