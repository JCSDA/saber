/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/generic/AtlasInterpolator.h"
#include "src/Covariance.h"
#include "src/Geometry.h"
#include "src/Increment.h"
#include "src/LinearVariableChange.h"
#include "src/ModelData.h"
#include "src/State.h"
#include "src/VariableChange.h"

namespace oops {
class AtlasInterpolator;
}  // namespace oops

namespace quench {

struct Traits {
  static std::string name()
    {return "quench";}
  static std::string nameCovar()
    {return "quenchCovariance";}

  typedef quench::Covariance           Covariance;
  typedef quench::Geometry             Geometry;
  typedef quench::Increment            Increment;
  typedef quench::LinearVariableChange LinearVariableChange;
  typedef quench::ModelData            ModelData;
  typedef quench::State                State;
  typedef quench::VariableChange       VariableChange;
};

struct Traits2 {
  static std::string name()
    {return "quench2";}
  static std::string nameCovar()
    {return "quenchCovariance2";}

  typedef quench::Covariance           Covariance;
  typedef quench::Geometry             Geometry;
  typedef quench::Increment            Increment;
  typedef quench::LinearVariableChange LinearVariableChange;
  typedef quench::ModelData            ModelData;
  typedef quench::State                State;
  typedef quench::VariableChange       VariableChange;
};

}  // namespace quench
