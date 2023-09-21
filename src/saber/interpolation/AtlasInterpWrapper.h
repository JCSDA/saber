/*
 * (C) Crown Copyright 2021-2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Point.h"
#include "atlas/util/PolygonLocator.h"
#include "atlas/util/PolygonXY.h"

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"


namespace atlas {
class Field;
class FieldSet;
}

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

class AtlasInterpWrapper {
 public:
  static const std::string classname() {return "saber::interpolation::AtlasInterpWrapper";}

  AtlasInterpWrapper(const atlas::grid::Partitioner &,
                     const atlas::FunctionSpace &,
                     const atlas::Grid &,
                     const atlas::FunctionSpace &);
  ~AtlasInterpWrapper() {}

  void execute(const atlas::Field &, atlas::Field &) const;

  void executeAdjoint(atlas::Field &, const atlas::Field &) const;

  void execute(const atlas::FieldSet &, atlas::FieldSet &) const;

  void executeAdjoint(atlas::FieldSet &, const atlas::FieldSet &) const;

 private:
  atlas::FunctionSpace localDstFunctionSpace_;
  std::vector<size_t> localTask_;
  atlas::FunctionSpace targetFunctionSpace_;
  atlas::Interpolation interp_;
  atlas::Redistribution redistr_;
  atlas::Redistribution inverseRedistr_;
};

}  // namespace interpolation
}  // namespace saber


