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
#include "oops/base/Variables.h"

#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierGridToSpectralParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierGridToSpectralParameters, SaberBlockParametersBase)

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierGridToSpectral : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierGridToSpectral";}

  typedef BifourierGridToSpectralParameters Parameters_;

  BifourierGridToSpectral(const oops::GeometryData &,
                          const oops::Variables &,
                          const eckit::Configuration &,
                          const Parameters_ &,
                          const oops::FieldSet3D &,
                          const oops::FieldSet3D &);
  virtual ~BifourierGridToSpectral() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return innerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override
    {}

 private:
  // Inner variables
  const oops::Variables innerVars_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransform> trans_;

  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
