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

#include "saber/bifourier/BifourierTransformBase.h"
#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/BlockParametersBase.h"
#include "saber/blocks/OuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierSpectralToGridParameters : public BlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierSpectralToGridParameters, BlockParametersBase)

 public:
  // Transform parameters
  oops::Parameter<BifourierTransformParameters> transform{"transform",
    BifourierTransformParameters(), this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierSpectralToGrid : public OuterBlockBase {
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
  void leftInverseMultiply(oops::FieldSet3D & fset) const override
    {inverseMultiply(fset);}
  void rightInverseMultiply(oops::FieldSet3D & fset) const override
    {inverseMultiply(fset);}

  void read() override
    {}

 private:
  // Inner geometry data
  std::unique_ptr<oops::GeometryData> innerGeometryData_;

  // Inner variables
  const oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  const std::shared_ptr<BifourierTransformBase> trans_;

  // Private methods

  // Inverse multiply
  void inverseMultiply(oops::FieldSet3D &) const;

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
