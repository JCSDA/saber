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
#include "oops/util/parameters/Parameters.h"

#include "saber/bifourier/BiperiodizationImpl.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BiperiodizationReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BiperiodizationReadParameters, oops::Parameters)

 public:
  // Input test file
  oops::RequiredParameter<eckit::LocalConfiguration> inputTestFile{"input test file", this};

  // Output inner test file
  oops::RequiredParameter<eckit::LocalConfiguration> outputInnerTestFile{"output inner test file",
    this};

  // Output outer test file
  oops::RequiredParameter<eckit::LocalConfiguration> outputOuterTestFile{"output outer test file",
    this};
};

// -----------------------------------------------------------------------------

class BiperiodizationParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BiperiodizationParameters, SaberBlockParametersBase)

 public:
  // Read parameters
  oops::OptionalParameter<BiperiodizationReadParameters> read{"read", this};

  // Biperiodization parameters
  oops::RequiredParameter<BiperiodizationImplParameters> biperParams{"biperiodization", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class Biperiodization : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::Biperiodization";}

  typedef BiperiodizationParameters Parameters_;

  Biperiodization(const oops::GeometryData &,
                  const oops::Variables &,
                  const eckit::Configuration &,
                  const Parameters_ &,
                  const oops::FieldSet3D &,
                  const oops::FieldSet3D &);
  virtual ~Biperiodization() = default;

  const oops::GeometryData & innerGeometryData() const override;
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  void read() override;

  void write() const override;

 private:
  // Inner grid
  atlas::StructuredGrid innerGrid_;

  // Inner partition
  std::vector<int> innerPartition_;

  // Inner FunctionSpace
  atlas::functionspace::StructuredColumns innerFs_;

  // Inner geometry data
  std::unique_ptr<oops::GeometryData> innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  const oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Biperiodization implementation
  std::unique_ptr<BiperiodizationImpl> biper_;

  // Test FieldSet3D
  std::unique_ptr<oops::FieldSet3D> inputTestFset_;
  std::unique_ptr<oops::FieldSet3D> outputInnerTestFset_;
  std::unique_ptr<oops::FieldSet3D> outputOuterTestFset_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
