/*
 * (C) Copyright 2022 UCAR
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

#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class EnsembleParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(EnsembleParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> localization{
        "localization", this};
  oops::Variables mandatoryActiveVars() const {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class Ensemble : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::Ensemble";}

  typedef EnsembleParameters Parameters_;

  Ensemble(const oops::GeometryData &,
           const std::vector<size_t> &,
           const oops::Variables &,
           const Parameters_ &,
           const atlas::FieldSet &,
           const atlas::FieldSet &,
           const std::vector<atlas::FieldSet> &);

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  std::vector<atlas::FieldSet> ensemble_;
  std::unique_ptr<SaberCentralBlockBase> loc_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
