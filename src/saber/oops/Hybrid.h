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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class CovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters, oops::Parameters)
 public:
  // Central and outer blocks
  oops::RequiredParameter<eckit::LocalConfiguration> saberCentralBlockConf{
    "saber central block", this};
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> saberOuterBlocksConf{
    "saber outer blocks", this};

  // Geometry [optional]
  oops::OptionalParameter<eckit::LocalConfiguration> hybridGeometry{"geometry", this};

  // Ensemble
  oops::Parameter<bool> iterativeEnsembleLoading{"iterative ensemble loading", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};
};

// -----------------------------------------------------------------------------

class WeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WeightParameters, oops::Parameters)
 public:
  // Scalar weight
  oops::OptionalParameter<double> value{"value", this};

  // File-base weight
  oops::OptionalParameter<eckit::LocalConfiguration> file{"file", this};
};

// -----------------------------------------------------------------------------

class ComponentParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ComponentParameters, oops::Parameters)
 public:
  // Covariance
  oops::RequiredParameter<CovarianceParameters> covariance{"covariance", this};

  // Weight
  oops::RequiredParameter<WeightParameters> weight{"weight", this};
};

// -----------------------------------------------------------------------------

class HybridParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(HybridParameters, SaberBlockParametersBase)
 public:
  // Vector of components
  oops::RequiredParameter<std::vector<ComponentParameters>> components{"components", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class Hybrid : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::Hybrid";}

  typedef HybridParameters Parameters_;

  Hybrid(const oops::GeometryData &,
         const oops::Variables &,
         const eckit::Configuration &,
         const Parameters_ & params,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const util::DateTime &,
         const size_t &)
      : SaberCentralBlockBase(params)
    {ABORT("the Hybrid block is a fake block, it should not be constructed");}

  virtual ~Hybrid() {}

  void randomize(atlas::FieldSet &) const override
    {ABORT("the Hybrid block is a fake block, it should not be used for randomization");}
  void multiply(atlas::FieldSet &) const override
    {ABORT("the Hybrid block is a fake block, it should not be used for multiplication");}

 private:
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
