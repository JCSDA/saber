/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

// Value or profile elemental parameters
class ValueOrProfileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ValueOrProfileParameters, oops::Parameters)

 public:
  // Group
  oops::RequiredParameter<std::string> group{"group", this};
  // Value
  oops::OptionalParameter<double> value{"value", this};
  // Profile
  oops::OptionalParameter<std::vector<double>> profile{"profile", this};
};

// -----------------------------------------------------------------------------

// Alias elemental paramaters
class AliasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AliasParameters, oops::Parameters)

 public:
  // In code
  oops::RequiredParameter<std::string> inCode{"in code", this};
  // In model file
  oops::RequiredParameter<std::string> inModelFile{"in model file", this};
};

// -----------------------------------------------------------------------------

class GroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GroupParameters, oops::Parameters)

 public:
  // Group name
  oops::RequiredParameter<std::string> name{"group name", this};

  // Variable in model file
  oops::OptionalParameter<std::string> varInModelFile{"variable in model file", this};

  // Group variables
  oops::RequiredParameter<std::vector<std::string>> variables{"variables", this};
};

// -----------------------------------------------------------------------------

class FastLAMParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FastLAMParametersBase, oops::Parameters)

 public:
  // Input files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> inputModelFilesConf{
    "input model files", this};

  // Horizontal length-scale
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>>
    rhFromYaml{"horizontal length-scale", this};

  // Vertical length-scale
  oops::OptionalParameter<std::vector<ValueOrProfileParameters>>
    rvFromYaml{"vertical length-scale", this};

  // Number of layers
  oops::OptionalParameter<size_t> nLayers{"number of layers", this};

  // Target resolution
  oops::OptionalParameter<size_t> resol{"resolution", this};

  // Skip tests
  oops::Parameter<bool> skipTests{"skip tests", false, this};

  // Interpolation accuracy tolerance
  oops::Parameter<double> accuracyTolerance{"interpolation accuracy tolerance", 0.05, this};

  // Interpolation adjoint tolerance
  oops::Parameter<double> adjointTolerance{"interpolation adjoint tolerance", 1.0e-12, this};

  // Data file (input / output)
  oops::OptionalParameter<std::string> dataFile{"data file", this};

  // Output model files
  oops::OptionalParameter<std::vector<eckit::LocalConfiguration>> outputModelFilesConf{
    "output model files", this};

  // Normalization accuracy stride (to reduce cost)
  oops::Parameter<size_t> normAccStride{"normalization accuracy stride", 10, this};

  // Multivariate strategy ('univariate', 'duplicated'), TODO(Benjamin): 'crossed'
  oops::RequiredParameter<std::string> strategy{"multivariate strategy", this};

  // Groups of variables
  oops::OptionalParameter<std::vector<GroupParameters>> groups{"groups", this};

  // Level for 2D variables ('first' or 'last')
  oops::Parameter<std::string> lev2d{"level for 2d variables", "first", this};

  // Aliases for model files
  oops::Parameter<std::vector<AliasParameters>> alias{"alias", {}, this};
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber

