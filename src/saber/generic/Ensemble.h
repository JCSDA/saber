/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class InflationFieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InflationFieldParameters, Parameters)

 public:
  // ATLAS inflation file
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFileConf{"atlas file", this};
  // Model inflation file
  oops::OptionalParameter<eckit::LocalConfiguration> modelFileConf{"model file", this};
};

// -----------------------------------------------------------------------------

class EnsembleParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(EnsembleParameters, SaberBlockParametersBase)

 public:
  // Inflation fields
  oops::OptionalParameter<InflationFieldParameters> inflationField{"inflation field", this};

  // Inflation value
  oops::Parameter<double> inflationValue{"inflation value", 1.0, this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------
// TODO(AS): remove this block once the parameters are refactored and passed
// correctly to SaberEnsembleBlockChain.
class Ensemble : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::Ensemble";}

  typedef EnsembleParameters Parameters_;

  Ensemble(const oops::GeometryData &,
           const oops::Variables &,
           const eckit::Configuration &,
           const Parameters_ & params,
           const oops::FieldSet3D & xb,
           const oops::FieldSet3D &) : SaberCentralBlockBase(params, xb.validTime())
  {throw eckit::Exception("the Ensemble block is a fake block, it should not be constructed",
    Here());}

  void randomize(oops::FieldSet3D &) const override
    {throw eckit::Exception("the Ensemble block is a fake block, it should not be used for"
      " randomization", Here());}
  void multiply(oops::FieldSet3D &) const override
    {throw eckit::Exception("the Ensemble block is a fake block, it should not be used for"
      " multiply", Here());}

 private:
  void print(std::ostream &) const override {}
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
