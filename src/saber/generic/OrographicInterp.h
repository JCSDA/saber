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

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class OrographicInterpParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(OrographicInterpParameters, SaberBlockParametersBase)

 public:
  // Orographic factor
  oops::Parameter<double> orographicFactor{"orographic factor", 1.0, this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class OrographicInterp : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::generic::OrographicInterp";}

  typedef OrographicInterpParameters Parameters_;

  OrographicInterp(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~OrographicInterp() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return innerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  const oops::Variables & innerVars_;

  // Parameters
  OrographicInterpParameters params_;

  // Interpolated variables
  oops::Variables interpVars_;

  // Interpolation FieldSet
  atlas::FieldSet interpFset_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
