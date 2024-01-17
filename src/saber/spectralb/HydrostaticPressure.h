/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/spectralb/GaussUVToGP.h"
#include "saber/vader/GpToHp.h"
#include "saber/vader/PressureParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {
  class GpToHp;
}

namespace spectralb {
  class GaussUVToGP;

// -----------------------------------------------------------------------------
/// \brief a saber block that creates hydrostatic pressure using
///        horizontal winds and unbalanced pressure. The "hydrostatic" part is
///        implied and not enforced in this saber block.
///        Note also that this saber block expects the outer and inner functionspaces
///        to be on the same Gaussian mesh.
///        To do this it uses GaussUVToGp and GpToHp saber blocks
class HydrostaticPressure : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::HydrostaticPressure";}

  typedef HydrostaticPressureParameters Parameters_;

  HydrostaticPressure(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~HydrostaticPressure();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  const oops::Variables innerVars_;
  const oops::Variables intermediateTempVars_;
  /// Gaussian (outer) functionspace
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  std::unique_ptr<saber::vader::GpToHp> gptohp_;
  std::unique_ptr<GaussUVToGP> gaussuvtogp_;
};

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
