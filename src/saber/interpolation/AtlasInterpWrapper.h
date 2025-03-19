/*
 * (C) Crown Copyright 2021-2022, Met Office
 * (C) Copyright 2024 Meteorologisk Institutt
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
#include "atlas/interpolation.h"
#include "atlas/redistribution/Redistribution.h"

#ifdef ATLAS_MAKE_SPARSE
#include "atlas/linalg/sparse.h"
#else
#include "eckit/linalg/SparseMatrix.h"
#endif

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

namespace atlas {
  class Field;
  class FieldSet;
  class Grid;
  namespace grid {
    class Partitioner;
  }
}

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

class AtlasInterpWrapper {
 public:
  static const std::string classname()
    {return "quench::interpolation::AtlasInterpWrapper";}

  AtlasInterpWrapper(const atlas::grid::Partitioner &,
                     const atlas::FunctionSpace &,
                     const atlas::Grid &,
                     const atlas::FunctionSpace &,
                     const std::string & interpType = "");
  ~AtlasInterpWrapper() {}

  void execute(const atlas::FieldSet &,
               atlas::FieldSet &) const;
  void executeAdjoint(atlas::FieldSet &,
                      const atlas::FieldSet &) const;

  eckit::linalg::SparseMatrix getInterpolationMatrix() const {
#ifdef ATLAS_MAKE_SPARSE
    return atlas::linalg::make_eckit_sparse_matrix(
                          atlas::interpolation::MatrixCache(interp_).matrix());
#else
    return atlas::interpolation::MatrixCache(interp_).matrix();
#endif
  }

  const atlas::FunctionSpace & getIntermediateFunctionSpace() const {
    return targetFspace_;
  }

  const atlas::Redistribution & getRedistribution() const {
    return redistr_;
  }

 private:
  void execute(const atlas::Field &,
               atlas::Field &) const;
  void executeAdjoint(atlas::Field &,
                      const atlas::Field &) const;

  atlas::FunctionSpace targetFspace_;
  atlas::Interpolation interp_;
  atlas::Redistribution redistr_;
  atlas::Redistribution inverseRedistr_;
};

}  // namespace interpolation
}  // namespace saber
