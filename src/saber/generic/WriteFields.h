/*
 * (C) Crown Copyright 2023 Met Office
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
namespace generic {

// -----------------------------------------------------------------------------

class WriteFieldsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(WriteFieldsParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}

  /// Path of output file.
  oops::RequiredParameter<std::string> outputPath{"output path", this};

  /// Label to assign to output file.
  oops::RequiredParameter<std::string> label{"label", this};

  /// Save fields to netCDF files.
  /// Whatever the value of this parameter, some basic information about each field is written
  /// to the test output stream.
  oops::Parameter<bool> saveNetCDFFile{"save netCDF file", true, this};

  /// List of fields to write out.
  /// If this parameter is empty (the default) then all fields contained in the
  /// FieldSet are written out.
  oops::Parameter<std::vector<std::string>> fieldNames{"field names", {}, this};

  /// Write out fields, contained in xb, in the block's constructor.
  oops::Parameter<bool> writeXb{"write xb", false, this};

  /// Write out fields, contained in fg, in the block's constructor.
  oops::Parameter<bool> writeFg{"write fg", false, this};

  /// Write out fields in the block's multiply() routine.
  oops::Parameter<bool> writeMultiply{"write multiply fset", false, this};

  /// Write out fields in the block's multiplyAD() routine.
  oops::Parameter<bool> writeMultiplyAD{"write multiplyad fset", false, this};

  /// Write out fields in the block's leftInverseMultiply() routine.
  oops::Parameter<bool> writeLeftInverseMultiply{"write leftinversemultiply fset", false, this};
};

// -----------------------------------------------------------------------------

class WriteFields : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::WriteFields";}

  typedef WriteFieldsParameters Parameters_;

  WriteFields(const oops::GeometryData &,
              const oops::Variables &,
              const eckit::Configuration &,
              const Parameters_ &,
              const oops::FieldSet3D &,
              const oops::FieldSet3D &);

  virtual ~WriteFields() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  /// Write selected fields to a file.
  void writeToFile(const oops::FieldSet3D &, const std::string &, size_t &) const;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const Parameters_ params_;
  mutable size_t count_xb_;
  mutable size_t count_fg_;
  mutable size_t count_multiply_;
  mutable size_t count_multiplyad_;
  mutable size_t count_leftinversemultiply_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
