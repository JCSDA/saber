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
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/bump/BUMPParameters.h"

#include "saber/bump/lib/BUMP.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class PsiChiToUVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(PsiChiToUVParameters, SaberBlockParametersBase)

 public:
  oops::OptionalParameter<BUMPParameters> readParams{"read", this};
  oops::OptionalParameter<BUMPParameters> calibrationParams{"calibration", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables({
    "stream_function",
    "velocity_potential",
    "eastward_wind",
    "northward_wind"});}
};

// -----------------------------------------------------------------------------

class PsiChiToUV : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::bump::PsiChiToUV";}

  typedef PsiChiToUVParameters Parameters_;

  PsiChiToUV(const oops::GeometryData &,
             const std::vector<size_t> &,
             const oops::Variables &,
             const eckit::Configuration &,
             const Parameters_ &,
             const atlas::FieldSet &,
             const atlas::FieldSet &);
  virtual ~PsiChiToUV();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;

  void read() override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  oops::Variables outerVars_;
  size_t levels_;
  BUMPParameters bumpParams_;
  std::unique_ptr<bump_lib::BUMP> bump_;
  size_t memberIndex_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
