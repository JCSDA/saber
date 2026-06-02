/*
 * (C) Crown Copyright 2026- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/interpolation/Interpolation.h"
#include "atlas/redistribution/Redistribution.h"

#pragma once

namespace saber {
namespace interpolation {

class SMVInterpWrapper {
 public:
  static const std::string classname() {return "quench::interpolation::SMVInterpWrapper";}

  SMVInterpWrapper(
      const atlas::functionspace::NodeColumns& CSFunctionSpace,
      const atlas::functionspace::StructuredColumns& gaussFunctionSpace,
      const atlas::Grid& gaussGrid,
      const bool includingVectorInterpolation = false,
      const bool normalizingInterpolation = false);
  ~SMVInterpWrapper() {}

  void doInterpolation(const std::vector<atlas::util::Config>& fieldConfigs,
                       const atlas::FieldSet& srcFieldSet,
                       atlas::FieldSet& newFieldSet) const;

 private:
  const atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  atlas::functionspace::NodeColumns gaussToCSFunctionSpace_;
  atlas::Interpolation interpolation_;
  atlas::Redistribution redistribution_;
  const bool includingVectorInterpolation_;
  const bool normalizingInterpolation_;
};

}  // namespace interpolation
}  // namespace saber
