/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/InterpolatorBase.h"
#include "saber/interpolation/InterpolatorBump.h"
#include "saber/interpolation/interpolatorbump_f.h"

namespace saber {

// This is the preferred place to define the static bump interpolator factory.
// However, strangely, it did not work with gcc (7 and 9.3) on ubuntu 18.04
// So, this has been moved to the header file.
//static oops::InterpolatorMaker<InterpolatorBump> makerBumpInterpolator_("bump");

// -----------------------------------------------------------------------------
InterpolatorBump::InterpolatorBump(const eckit::Configuration & config,
                   const atlas::FunctionSpace & infspace,
                   const atlas::FunctionSpace & outfspace,
                   const atlas::field::FieldSetImpl * masks,
                   const eckit::mpi::Comm & comm) {
  bint_create_f90(keyBumpInterpolator_, &comm, infspace.get(), outfspace.get(),
                  masks, config);
}

// -----------------------------------------------------------------------------
  void InterpolatorBump::apply(atlas::FieldSet const & infields,
                               atlas::FieldSet & outfields) {
  bint_apply_f90(keyBumpInterpolator_, infields.get(), outfields.get());
}

// -----------------------------------------------------------------------------
  void InterpolatorBump::apply_ad(atlas::FieldSet const & fields_grid2,
                                  atlas::FieldSet & fields_grid1) {
  bint_apply_ad_f90(keyBumpInterpolator_, fields_grid2.get(),
                    fields_grid1.get());
}

// -----------------------------------------------------------------------------
InterpolatorBump::~InterpolatorBump() {
  bint_delete_f90(keyBumpInterpolator_);
}

// -----------------------------------------------------------------------------
void InterpolatorBump::print(std::ostream & os) const {
  os << " InterpolatorBump: print not implemented yet.";
}

// -----------------------------------------------------------------------------
}  // namespace saber
