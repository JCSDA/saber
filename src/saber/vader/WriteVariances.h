/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/util/Calibration.h"
#include "saber/vader/VarianceAccumulationUtils.h"

namespace saber {
namespace vader {

class calibrationVariancesParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(calibrationVariancesParameters, oops::Parameters)

 public:
  oops::RequiredParameter<util::calibrationWriteParameters> writeParams{"write", this};
};

class binningParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(binningParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> type{"type", this};

  // optional parameters needed for writing file
  oops::OptionalParameter<std::string> mpiRankPattern{"mpi rank pattern", this};

  oops::Parameter<bool> oneFilePerTask{"one file per task", true, this};

  oops::OptionalParameter<std::string> filePath{"file path", this};

  // optional parameters needed for latitude band binning
  oops::OptionalParameter<std::size_t> noOfBins{"no of bins", this};

  oops::OptionalParameter<std::vector<double>> lowerBounds{"lower bounds", this};

  oops::OptionalParameter<std::vector<double>> upperBounds{"upper bounds", this};

  oops::Parameter<bool> includeLowerBound{"include lower bound", false, this};

  oops::Parameter<bool> includeUpperBound{"include upper bound", true, this};

  oops::Parameter<bool> includeOuterBound{"include outer bound", true, this};

  oops::OptionalParameter<std::string> stateFieldName{"state field name", this};
};

class varianceInstantaneousParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(varianceInstantaneousParameters, oops::Parameters)

 public:
  /// Write out fields in the block's multiply() routine with filename below
  oops::OptionalParameter<std::string> multiplyFileName{"multiply fset filename", this};

  /// Write out fields in the block's multiplyAD() routine with filename below
  oops::OptionalParameter<std::string> multiplyADFileName{"multiplyad fset filename", this};

  /// Write out fields in the block's leftInverseMultiply() routine with filename below
  oops::OptionalParameter<std::string> leftInverseFileName{"left inverse fset filename", this};

  oops::RequiredParameter<std::string> outputPath{"output path", this};
};

class WriteVariancesParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(WriteVariancesParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}

  oops::OptionalParameter<calibrationVariancesParameters>
    calibrationParams{"calibration", this};

  /// Save fields to netCDF files.
  /// Whatever the value of this parameter, some basic information about each field is written
  /// to the test output stream.
  oops::Parameter<bool> saveNetCDFFile{"save netCDF file", true, this};

  oops::Parameter<std::string> statisticsType{"statistics type",
                                              "variance", this};

  /// List of field names to specify fields to write to file
  /// When used with instantanous statistics - if empty all Fields contained in fieldset
  /// considered.
  /// When used in calibration mode - lists all fields that involve variances
  /// / vertical covariances with same field names ... inter-variable cross covariances are
  /// separately determined.
  oops::Parameter<std::vector<std::string>> fieldNames{"field names", {}, this};

  oops::RequiredParameter<binningParameters> binning{"binning", this};

  oops::OptionalParameter<varianceInstantaneousParameters>
    instantaneousParams{"instantaneous statistics", this};

  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>>
    crosscovParams{"additional cross covariances", this};
};

// -----------------------------------------------------------------------------
// Calculates the global-averaged variances for each model level and field,
// It is to be used mainly as a diagnostic saber block.
// It currently works for regular Gaussian and cubed-sphere dual grids.
// In the future it will be extended to latitude bands and maybe grid-point variances


class WriteVariances : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::WriteVariances";}

  typedef WriteVariancesParameters Parameters_;

  WriteVariances(const oops::GeometryData &,
                 const oops::Variables &,
                 const eckit::Configuration &,
                 const Parameters_ &,
                 const oops::FieldSet3D &,
                 const oops::FieldSet3D &);

  virtual ~WriteVariances() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}

  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;
  void directCalibration(const oops::FieldSets &) override;
  void write() const override;

 private:
  void print(std::ostream &) const override;

  void writeInstantVariances(const eckit::mpi::Comm & comm,
                             const atlas::FieldSet & fset,
                             const std::string & description) const;

  void diagnostics(const std::string & tag,
                   const oops::FieldSet3D & fset) const;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const Parameters_ params_;
  std::string binType_;
  eckit::LocalConfiguration netCDFConf_;
  std::size_t sizeOwned_;
  std::size_t totalBins_;
  atlas::FieldSet binningData_;   // holds the latitude, longitude and
                                  // binning area-weights and index information
  atlas::FieldSet ensembleStats_;  // in calibration mode
                                   // has the statistics
  const util::DateTime datetime_;
  mutable std::size_t count_multiply_;
  mutable std::size_t count_multiplyad_;
  mutable std::size_t count_leftinversemultiply_;
  std::size_t sample_size_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
