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

class BifourierSpectralToGridParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierSpectralToGridParameters, SaberBlockParametersBase)

 public:
  // Truncation type
  oops::Parameter<std::string> truncationType{"truncation type", "arome", this};

  // Skip tests
  oops::Parameter<bool> skipTests{"skip tests", false, this};

  // Spectral tests tolerance
  oops::Parameter<double> specTolerance{"spectral tolerance", 1.0e-9, this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierSpectralToGrid : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierSpectralToGrid";}

  typedef BifourierSpectralToGridParameters Parameters_;

  BifourierSpectralToGrid(const oops::GeometryData &,
                          const oops::Variables &,
                          const eckit::Configuration &,
                          const Parameters_ &,
                          const oops::FieldSet3D &,
                          const oops::FieldSet3D &);
  virtual ~BifourierSpectralToGrid() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return *innerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override
    {}

 private:
  // Inner geometry data
  std::unique_ptr<oops::GeometryData> innerGeometryData_;

  // Inner variables
  const oops::Variables innerVars_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransform> trans_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
