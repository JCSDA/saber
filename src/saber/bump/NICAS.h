/*
 * (C) Copyright 2021 UCAR
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

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class NICASParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(NICASParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMPParameters> bumpParams{"bump", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class NICAS : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::bump::NICAS";}

  typedef NICASParameters Parameters_;

  NICAS(const oops::GeometryData &,
        const std::vector<size_t> &,
        const oops::Variables &,
        const Parameters_ &,
        const atlas::FieldSet &,
        const atlas::FieldSet &,
        const std::vector<atlas::FieldSet> &,
        const size_t &);
  virtual ~NICAS();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
