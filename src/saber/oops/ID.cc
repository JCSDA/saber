/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/ID.h"

#include <memory>
#include <string>
#include <vector>

#include "oops/util/Timer.h"

#include "oops/base/GeometryData.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<ID> makerID_("ID");

// -----------------------------------------------------------------------------

ID::ID(const eckit::mpi::Comm & comm,
       const oops::GeometryData & geometryData,
       const std::vector<size_t> & activeVariableSizes,
       const eckit::Configuration & conf,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg,
       const std::vector<atlas::FieldSet> & fsetVec)
  : SaberCentralBlockBase(conf)
{
  oops::Log::trace() << classname() << "::ID starting" << std::endl;

  // Deserialize configuration
  IDParameters params;
  params.validateAndDeserialize(conf);

  oops::Log::trace() << classname() << "::ID done" << std::endl;
}

// -----------------------------------------------------------------------------

ID::~ID() {
  oops::Log::trace() << classname() << "::~ID starting" << std::endl;
  util::Timer timer(classname(), "~ID");
  oops::Log::trace() << classname() << "::~ID done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber
