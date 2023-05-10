/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/Spectral.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

class SpectralToSpectralParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SpectralToSpectralParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<int> inputTruncation{"input truncation", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};


class SpectralToSpectral : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SpectralToSpectral";}

  typedef SpectralToSpectralParameters Parameters_;

  SpectralToSpectral(const oops::GeometryData &,
                     const std::vector<size_t> &,
                     const oops::Variables &,
                     const eckit::Configuration &,
                     const SpectralToSpectralParameters &,
                     const atlas::FieldSet &,
                     const atlas::FieldSet &,
                     const util::DateTime &);

  virtual ~SpectralToSpectral() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return activeVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  void truncate_or_extend(const atlas::functionspace::Spectral &,
                          const atlas::functionspace::Spectral &,
                          atlas::FieldSet &) const;
  void truncate(const atlas::functionspace::Spectral &,
                const atlas::functionspace::Spectral &,
                atlas::FieldSet &) const;
  void truncateAD(const atlas::functionspace::Spectral &,
                  const atlas::functionspace::Spectral &,
                  atlas::FieldSet &) const;

  const atlas::functionspace::Spectral innerFunctionSpace_;
  const atlas::functionspace::Spectral outerFunctionSpace_;
  const oops::GeometryData innerGeometryData_;
  const oops::Variables activeVars_;
};

}  // namespace spectralb
}  // namespace saber
