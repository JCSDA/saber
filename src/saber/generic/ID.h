/*
 * (C) Copyright 2021 UCAR
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
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class IDParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(IDParameters, SaberBlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class IDCentral : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::IDCentral";}

  typedef IDParameters Parameters_;

  IDCentral(const oops::GeometryData &,
            const oops::Variables &,
            const eckit::Configuration &,
            const Parameters_ &,
            const oops::FieldSet3D &,
            const oops::FieldSet3D &);

  virtual ~IDCentral();

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  size_t ctlVecSize() const override {return ctlVecSize_;}
  void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const override;

 private:
  const oops::GeometryData & geometryData_;
  const oops::Variables activeVars_;
  size_t ctlVecSize_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

class IDOuter : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::generic::IDOuter";}

  typedef IDParameters Parameters_;

  IDOuter(const oops::GeometryData &,
          const oops::Variables &,
          const eckit::Configuration &,
          const Parameters_ &,
          const oops::FieldSet3D &,
          const oops::FieldSet3D &);

  virtual ~IDOuter() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
