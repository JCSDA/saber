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
  } else {
    throw eckit::Exception("wrong interpolation type", Here());
  }

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::execute(const atlas::FieldSet & srcFieldSet,
                            atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->execute(srcFieldSet, targetFieldSet);
  }

  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                   const atlas::FieldSet & targetFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->executeAdjoint(srcFieldSet, targetFieldSet);
  }

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::insertVerticalInterpolation(const std::string & var,
                                                const std::vector<InterpElement> & item) {
  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;

  if (verInterps_.find(var) != verInterps_.end()) {
    throw eckit::Exception("vertical interpolation already computed for this variables");
  }
  verInterps_.insert({var, item});

  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
