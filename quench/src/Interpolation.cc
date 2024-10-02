/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Interpolation.h"

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"

// -----------------------------------------------------------------------------

namespace quench {

// -----------------------------------------------------------------------------

Interpolation::Interpolation(const eckit::Configuration & conf,
                             const eckit::mpi::Comm & comm,
                             const atlas::grid::Partitioner & srcPartitioner,
                             const atlas::FunctionSpace & srcFspace,
                             const std::string & srcUid,
                             const atlas::Grid & dstGrid,
                             const atlas::FunctionSpace & dstFspace,
                             const std::string & dstUid)
  : srcUid_(srcUid), dstUid_(dstUid), dstFspace_(dstFspace), atlasInterpWrapper_() {
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;

  // Get interpolation type
  const std::string type = conf.getString("interpolation type");

  // Setup interpolation
  if (type == "atlas interpolation wrapper") {
    atlasInterpWrapper_ = std::make_shared<saber::interpolation::AtlasInterpWrapper>(srcPartitioner,
      srcFspace, dstGrid, dstFspace);
  } else if (type == "regional") {
    regionalInterp_ = std::make_shared<atlas::Interpolation>(
      atlas::util::Config("type", "regional-linear-2d"),
      srcFspace, dstFspace);
  } else {
    throw eckit::Exception("wrong interpolation type", Here());
  }

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::execute(const atlas::FieldSet & srcFieldSet,
                            atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->execute(srcFieldSet, tgtFieldSet);
  }
  if (regionalInterp_) {
    regionalInterp_->execute(srcFieldSet, tgtFieldSet);
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                   const atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->executeAdjoint(srcFieldSet, tgtFieldSet);
  }
  if (regionalInterp_) {
    regionalInterp_->execute_adjoint(srcFieldSet, tgtFieldSet);
  }

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------


void Interpolation::insertVerticalInterpolation(const std::string & var,
                                                const std::vector<std::array<size_t, 2>> & stencil,
                                                const std::vector<std::array<double, 2>> & weights,
                                                const std::vector<size_t> & stencilSize) {
  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;

  if (verStencil_.find(var) != verStencil_.end()) {
    throw eckit::Exception("vertical interpolation already computed for this variables");
  }
  ASSERT(stencil.size() == stencilSize.size());
  ASSERT(weights.size() == stencilSize.size());
  verStencil_.insert({var, stencil});
  verWeights_.insert({var, weights});
  verStencilSize_.insert({var, stencilSize});

  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeVertical(const atlas::FieldSet & srcFieldSet,
                                    atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeVertical starting" << std::endl;

  for (auto & tgtField : tgtFieldSet) {
    const std::string var = tgtField.name();
    const auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[var]);
    auto tgtView = atlas::array::make_view<double, 1>(tgtField);
    tgtView.assign(0.0);
    for (size_t jo = 0; jo < verStencilSize_.at(var).size(); ++jo) {
      for (size_t jj = 0; jj < verStencilSize_.at(var)[jo]; ++jj) {
        tgtView(jo) += verWeights_.at(var)[jo][jj]*srcView(jo, verStencil_.at(var)[jo][jj]);
      }
    }
  }

  oops::Log::trace() << classname() << "::executeVertical done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeVerticalAdjoint(atlas::FieldSet & srcFieldSet,
                                           const atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeVerticalAdjoint starting" << std::endl;

  for (const auto & tgtField : tgtFieldSet) {
    const std::string var = tgtField.name();
    auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[var]);
    const auto tgtView = atlas::array::make_view<double, 1>(tgtField);
    srcView.assign(0.0);
    for (size_t jo = 0; jo < verStencilSize_.at(var).size(); ++jo) {
      for (size_t jj = 0; jj < verStencilSize_.at(var)[jo]; ++jj) {
        srcView(jo, verStencil_.at(var)[jo][jj]) += verWeights_.at(var)[jo][jj]*tgtView(jo);
      }
    }
  }

  oops::Log::trace() << classname() << "::executeVerticalAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
