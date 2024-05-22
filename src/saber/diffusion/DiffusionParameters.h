/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/util/parameters/ParameterTraits.h"
#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {

class DiffusionParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(DiffusionParameters, saber::SaberBlockParametersBase)

 public:
  // ----------------------------------------------------------------------------
  enum class Strategy {UNIVARIATE, DUPLICATED};

  struct StrategyParameterTraitsHelper {
    typedef Strategy EnumType;
    static constexpr char enumTypeName[] = "Strategy";
    static constexpr util::NamedEnumerator<Strategy> namedValues[] = {
      { Strategy::UNIVARIATE, "univariate"},
      { Strategy::DUPLICATED, "duplicated"}
    };
  };

  // ----------------------------------------------------------------------------

  /// parameters for a group of variables that will have the same
  /// correlation/localization applied.
  /// TODO(Travis) implement proper parameters for "horizontal", "vertical", and "write"
  class Group : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Group, oops::Parameters)
   public:
    oops::OptionalParameter<std::vector<std::string> > variables{"variables", this};
    oops::OptionalParameter<eckit::LocalConfiguration> horizontal{"horizontal", this};
    oops::OptionalParameter<eckit::LocalConfiguration> vertical{"vertical", this};
    oops::OptionalParameter<eckit::LocalConfiguration> write{"write", this};
    oops::Parameter<Strategy> multivariateStrategy{"multivariate strategy",
                                                   Strategy::UNIVARIATE, this};
  };

  // ----------------------------------------------------------------------------

  class Calibration : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Calibration, oops::Parameters)
   public:
    oops::RequiredParameter<int> normalizationIterations {"normalization.iterations", this};
    oops::RequiredParameter<std::vector<Group> > groups {"groups", this};
  };

  // ----------------------------------------------------------------------------

  class Read : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Read, oops::Parameters)
   public:
    oops::RequiredParameter<std::vector<Group> > groups {"groups", this};
  };

  // ----------------------------------------------------------------------------

  oops::OptionalParameter<Read> read{"read", this};
  oops::OptionalParameter<Calibration> calibration{"calibration", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};
}  // namespace saber

namespace oops {
  template <>
  struct ParameterTraits<saber::DiffusionParameters::Strategy> :
    public EnumParameterTraits<saber::DiffusionParameters::StrategyParameterTraitsHelper>{};
}

