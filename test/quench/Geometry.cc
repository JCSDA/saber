/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Geometry.h"

#include <math.h>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/projection.h"
#include "atlas/util/Point.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Geometry::Geometry(const Parameters_ & params,
                   const eckit::mpi::Comm & comm) : comm_(comm), levels_(1) {
  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  if (params.functionSpace.value() == "StructuredColumns") {
    // StructuredColumns function space

    // Setup grid
    gridConfig_ = params.grid.value();
    oops::Log::info() << "Grid config: " << gridConfig_ << std::endl;
    atlasGrid_.reset(new atlas::Grid(gridConfig_));

    // Setup partitioner
    const atlas::grid::Partitioner partitioner(params.partitioner.value());

    // Setup distribution
    const atlas::grid::Distribution distribution(*atlasGrid_, partitioner);

    // Number of levels
    levels_ = params.levels.value();

    // Vertical unit
    const boost::optional<std::vector<double>> &vunit = params.vunit.value();
    for (size_t jlevel = 0; jlevel < levels_; ++jlevel) {
      if (vunit != boost::none) {
        vunit_.push_back((*vunit)[jlevel]);
      } else {
        vunit_.push_back(static_cast<double>(jlevel+1));
      }
    }

    // Setup function space
    atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(*atlasGrid_,
                              distribution));

    // Print summary
    this->print(oops::Log::info());
  } else {
    ABORT(params.functionSpace.value() + " function space not implemented yet");
  }

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  atlas::Field vunit = atlasFunctionSpace_->createField<double>(
    atlas::option::name("vunit") | atlas::option::levels(levels_));
  auto view = atlas::array::make_view<double, 2>(vunit);
  for (atlas::idx_t jnode = 0; jnode < vunit.shape(0); ++jnode) {
    for (atlas::idx_t jlevel = 0; jlevel < vunit.shape(1); ++jlevel) {
       view(jnode, jlevel) = vunit_[jlevel];
    }
  }
  atlasFieldSet_->add(vunit);
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_), levels_(other.levels_),
  vunit_(other.vunit_) {
  // Copy ATLAS grid
  gridConfig_ = other.gridConfig_;
  atlasGrid_.reset(new atlas::Grid(gridConfig_));

  // Copy ATLAS function space
  if (other.atlasFunctionSpace_->type() == "StructuredColumns") {
     atlasFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(*(
                               other.atlasFunctionSpace_)));
  } else {
    ABORT(other.atlasFunctionSpace_->type() + " function space not implemented yet");
  }

  // Copy ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  atlas::Field vunit = other.atlasFieldSet()->field("vunit");
  atlasFieldSet_->add(vunit);
}
// -------------------------------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes(vars.size(), levels_);
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  os << "Quench geometry grid:" << std::endl;
  os << "- name: " << atlasGrid_->name() << std::endl;
  os << "- size: " << atlasGrid_->size() << std::endl;
  os << "Function space:" << std::endl;
  os << "- type: " << atlasFunctionSpace_->type() << std::endl;
  os << "- size: " << atlasFunctionSpace_->size() << std::endl;
  os << "Vertical levels: " << std::endl;
  os << "- number: " << levels_ << std::endl;
  os << "- vunit: " << vunit_ << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace quench
