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
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierSpectralConverterParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierSpectralConverterParameters, SaberBlockParametersBase)

 public:
  // Background variable used to define the inner grid-poind function space
  oops::OptionalParameter<std::string> fspaceFromBkgVar{
    "inner grid-point function space from background variable", this};

  // Number of grid-points in X direction
  oops::OptionalParameter<int> nx{"nx", this};

  // Number of grid-points in Y direction
  oops::OptionalParameter<int> ny{"ny", this};

  // Partitioner
  oops::Parameter<std::string> partitioner{"partitioner", "checkerboard", this};

  // Transform parameters
  oops::Parameter<BifourierTransformParameters> transform{"transform",
    BifourierTransformParameters(), this};

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierSpectralConverter : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierSpectralConverter";}

  typedef BifourierSpectralConverterParameters Parameters_;

  BifourierSpectralConverter(const oops::GeometryData &,
                             const oops::Variables &,
                             const eckit::Configuration &,
                             const Parameters_ &,
                             const oops::FieldSet3D &,
                             const oops::FieldSet3D &);
  virtual ~BifourierSpectralConverter() = default;

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
  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Inner variables
  const oops::Variables innerVars_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  std::shared_ptr<BifourierTransformBase> innerTrans_;
  std::shared_ptr<BifourierTransformBase> outerTrans_;

  // Inner grid-point geometry data
  std::unique_ptr<oops::GeometryData> innerGpGeometryData_;

  // Inner spectral geometry data
  std::unique_ptr<oops::GeometryData> innerGeometryData_;

  // Communication
  size_t sendSize_;
  size_t recvSize_;
  std::vector<int> sendCounts_;
  std::vector<int> recvCounts_;
  std::vector<int> sendDispls_;
  std::vector<int> recvDispls_;
  std::vector<int> sendIndex_;
  std::vector<int> recvIndex_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
