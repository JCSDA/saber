/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_
#define SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/vader/HydrostaticExnerParameters.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
class HydrostaticExnerSaberBlockParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HydrostaticExnerSaberBlockParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<hydrostaticexnerParameters>
    hydrostaticexnerParams{"covariance data", this};
};

// -----------------------------------------------------------------------------
// This saber block is here to do 3 jobs
// 1) the vertical regression on geostrophic pressure
// 2) summing the result with unbalanced pressure to create hydrostatic_pressure
// 3) converting hydrostatic pressure to exner pressure.
// -----------------------------------------------------------------------------
class HydrostaticExnerSaberBlock : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::HydrostaticExnerSaberBlock";}

  typedef HydrostaticExnerSaberBlockParameters Parameters_;

  HydrostaticExnerSaberBlock(const atlas::FunctionSpace &,
                             const atlas::FieldSet &,
                             const std::vector<size_t> &,
                             const Parameters_ &,
                             const atlas::FieldSet &,
                             const atlas::FieldSet &,
                             const std::vector<atlas::FieldSet> &);
  virtual ~HydrostaticExnerSaberBlock();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  oops::Variables inputVars_;
  atlas::FieldSet covFieldSet_;
  atlas::FieldSet augmentedStateFieldSet_;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_VADER_HYDROSTATICEXNERSABERBLOCK_H_
