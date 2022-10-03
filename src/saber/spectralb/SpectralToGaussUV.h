/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SpectralToGaussUVParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralToGaussUVParameters, SaberOuterBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> gaussGridUid{"gauss grid uid",
    "Gauss Grid UID", this};
  oops::Parameter<bool> useStreamFunctionVelocityPotential{
    "use streamfunction and velocity potential", false, this};
};

// -----------------------------------------------------------------------------


class SpectralToGaussUV : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralToGaussUV";}

  typedef SpectralToGaussUVParameters Parameters_;

  SpectralToGaussUV(const oops::GeometryData &,
         const std::vector<std::size_t> &,
         const oops::Variables &,
         const Parameters_ &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);

  virtual ~SpectralToGaussUV() = default;

  const oops::GeometryData & inputGeometryData() const override {return inputGeometryData_;}
  const oops::Variables & inputVars() const override {return inputVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;

  Parameters_ params_;
  oops::Variables inputVars_;
  oops::Variables outputVars_;
  std::vector<std::size_t> activeVariableSizes_;

  /// Gaussian (output) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  /// Spectral (input) functionspace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Trans object for gaussian-spectral transforms
  const atlas::trans::Trans trans_;
  const oops::GeometryData inputGeometryData_;
};

}  // namespace spectralb
}  // namespace saber
