/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierSplitTPsParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierSplitTPsParameters, SaberBlockParametersBase)

 public:
  // Backward mode
  oops::Parameter<bool> backward{"backward mode", false, this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables(
    std::vector<std::string>({
    "air_temperature",
    "log_of_air_pressure_at_surface",
    "air_temperature_and_log_of_air_pressure_at_surface"}));}
};

// -----------------------------------------------------------------------------

class BifourierSplitTPs : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierSplitTPs";}

  typedef BifourierSplitTPsParameters Parameters_;

  BifourierSplitTPs(const oops::GeometryData &,
                    const oops::Variables &,
                    const eckit::Configuration &,
                    const Parameters_ &,
                    const oops::FieldSet3D &,
                    const oops::FieldSet3D &);
  virtual ~BifourierSplitTPs() = default;

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
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransform> trans_;

  // Number of levels
  size_t nz_;

  // Private methods

  // Forward application
  void forward(oops::FieldSet3D &) const;

  // Backward application
  void backward(oops::FieldSet3D &) const;

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
