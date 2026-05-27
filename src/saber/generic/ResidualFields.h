/*
 * (C) Crown Copyright 2026 Met Office
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
#include "oops/util/ParallelFieldSetIO.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class ResidualFieldsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ResidualFieldsParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}

  /// Path of input file.
  oops::RequiredParameter<std::string> inputPath{"input path", this};

  oops::Parameter<bool> readParallelIONetCDFFile{"parallel IO", false, this};

  /// List of fields to write out.
  /// If this parameter is empty (the default) then all fields contained in the
  /// FieldSet are written out.
  oops::Parameter<std::vector<std::string>> fieldNames{"field names", {}, this};

  /// Write out fields in the block's multiply() routine with filename below
  oops::Parameter<std::string> multiplyFileName{"multiply fset filename", {}, this};
};

// -----------------------------------------------------------------------------

class ResidualFields : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::ResidualFields";}

  typedef ResidualFieldsParameters Parameters_;

  ResidualFields(const oops::GeometryData &,
                 const oops::Variables &,
                 const eckit::Configuration &,
                 const Parameters_ &,
                 const oops::FieldSet3D &,
                 const oops::FieldSet3D &);

  virtual ~ResidualFields() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  const Parameters_ params_;
  std::unique_ptr<const util::ParallelFieldSetIO> io_;
  mutable size_t count_multiply_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
