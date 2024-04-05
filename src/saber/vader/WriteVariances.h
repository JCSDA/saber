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

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

class WriteVariancesParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(WriteVariancesParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}

  /// Path of output file.
  oops::RequiredParameter<std::string> outputPath{"output path", this};

  /// Save fields to netCDF files.
  /// Whatever the value of this parameter, some basic information about each field is written
  /// to the test output stream.
  oops::Parameter<bool> saveNetCDFFile{"save netCDF file", true, this};

  /// List of fields to write out.
  /// If this parameter is empty (the default) then all fields contained in the
  /// FieldSet are written out.
  oops::Parameter<std::vector<std::string>> fieldNames{"field names", {}, this};

  /// Write out fields in the block's multiply() routine with filename below
  oops::OptionalParameter<std::string> multiplyFileName{"multiply fset filename", this};

  /// Write out fields in the block's multiplyAD() routine with filename below
  oops::OptionalParameter<std::string> multiplyADFileName{"multiplyad fset filename", this};

  /// Write out fields in the block's leftInverseMultiply() routine with filename below
  oops::OptionalParameter<std::string> leftInverseFileName{"left inverse fset filename", this};
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

 private:
  void print(std::ostream &) const override;

  void writeToFile(const eckit::mpi::Comm & comm,
                   const atlas::FieldSet & fset,
                   const std::string & description) const;

  void diagnostics(const std::string & tag,
                   const oops::FieldSet3D & fset) const;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const Parameters_ params_;
  mutable size_t count_multiply_;
  mutable size_t count_multiplyad_;
  mutable size_t count_leftinversemultiply_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
