/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalAtlasInterpolator.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/interpolation/Geometry.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------
class InterpolationParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<GeometryParameters> innerGeom{"inner geometry", this};
  oops::Parameter<eckit::LocalConfiguration> localInterpConf{"local interpolator",
    eckit::LocalConfiguration(), this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class Interpolation : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const oops::GeometryData &,
                const oops::Variables &,
                const eckit::Configuration &,
                const Parameters_ &,
                const atlas::FieldSet &,
                const atlas::FieldSet &,
                const util::DateTime &);
  virtual ~Interpolation() = default;

  const oops::GeometryData & innerGeometryData() const override {return *innerGeomData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

  atlas::FieldSet generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                        const oops::Variables & innerVars,
                                        const size_t & timeRank) const override
    {return util::createSmoothFieldSet(innerGeometryData.comm(),
                                       innerGeometryData.functionSpace(),
                                       innerVars);}

  atlas::FieldSet generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                        const oops::Variables & outerVars,
                                        const size_t & timeRank) const override
    {return util::createSmoothFieldSet(outerGeometryData.comm(),
                                       outerGeometryData.functionSpace(),
                                       outerVars);}

 private:
  void print(std::ostream &) const override;

  const Parameters_ params_;
  const oops::GeometryData & outerGeomData_;
  const oops::Variables innerVars_;
  const oops::Variables activeVars_;
  // pointers for delayed initialization
  std::unique_ptr<oops::GeometryData> innerGeomData_;
  std::unique_ptr<oops::GlobalAtlasInterpolator> interp_;
  mutable std::unique_ptr<oops::GlobalAtlasInterpolator> inverseInterp_;
};

}  // namespace interpolation
}  // namespace saber
