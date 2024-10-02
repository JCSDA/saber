/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace quench {

class Covariance;
class Geometry;
class Increment;
class LinearVariableChange;
class ModelData;
class State;
class VariableChange;

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

}  // namespace quench
