/*
 * (C) Crown Copyright 2022-2023 Met Office
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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/PressureParameters.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------
/// \brief saber block that converts zonal and meridional wind to geostrophic
///        pressure on a Gaussian latitude mesh.

class GaussUVToGP : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::GaussUVToGP";}

  typedef GaussUVToGPParameters Parameters_;

  GaussUVToGP(const oops::GeometryData &,
              const oops::Variables &,
              const eckit::Configuration &,
              const Parameters_ &,
              const oops::FieldSet3D &,
              const oops::FieldSet3D &);

  virtual ~GaussUVToGP() = default;

  const oops::GeometryData & innerGeometryData()
    const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  Parameters_ params_;
  const oops::Variables outerVars_;
  const oops::Variables innerVars_;
  const oops::Variables activeOuterVars_;
  const oops::Variables innerOnlyVars_;

  /// Gaussian (outer) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  /// Spectral (inner) functionspace
  const atlas::functionspace::Spectral specFunctionSpace_;
  /// Trans object for gaussian-spectral transforms
  const atlas::trans::Trans trans_;
  const oops::GeometryData innerGeometryData_;
  const atlas::FieldSet augmentedState_;
};

}  // namespace spectralb
}  // namespace saber
