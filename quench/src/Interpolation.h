/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/interpolation.h"

#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "saber/interpolation/AtlasInterpWrapper.h"

namespace atlas {
  class Field;
  class Grid;
  namespace grid {
    class Partitioner;
  }
}

namespace quench {
  class Geometry;

// -----------------------------------------------------------------------------
/// Interpolation class

class Interpolation : public util::Printable,
                      private util::ObjectCounter<Interpolation> {
 public:
  static const std::string classname()
    {return "quench::Interpolation";}

  // Constructor/destructor
  Interpolation(const Geometry &,
                const Geometry &);
  ~Interpolation()
    {}

  // Horizontal interpolation and adjoint
  void execute(const atlas::FieldSet &,
               atlas::FieldSet &) const;
  void executeAdjoint(atlas::FieldSet &,
                      const atlas::FieldSet &) const;

  // Vertical interpolation
  void insertVerticalInterpolation(const std::string &,
                                   const std::vector<std::array<size_t, 2>> &,
                                   const std::vector<std::array<double, 2>> &,
                                   const std::vector<size_t> &);
  void executeVertical(const atlas::FieldSet &,
                       atlas::FieldSet &) const;
  void executeVerticalAdjoint(atlas::FieldSet &,
                              const atlas::FieldSet &) const;

 private:
  // Print
  void print(std::ostream &) const;

  // Destination function space
  atlas::FunctionSpace tgtFspace_;

  // ATLAS interpolation wrapper from SABER
  std::shared_ptr<saber::interpolation::AtlasInterpWrapper> atlasInterpWrapper_;

  // Regional ATLAS interpolation
  std::shared_ptr<atlas::Interpolation> regionalInterp_;

  // OOPS unstructured interpolation
  std::shared_ptr<oops::UnstructuredInterpolator> unstructuredInterp_;

  // Vertical interpolations
  std::unordered_map<std::string, std::vector<std::array<size_t, 2>>> verStencil_;
  std::unordered_map<std::string, std::vector<std::array<double, 2>>> verWeights_;
  std::unordered_map<std::string, std::vector<size_t>> verStencilSize_;
};

}  // namespace quench
