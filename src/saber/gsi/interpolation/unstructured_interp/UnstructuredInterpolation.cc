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
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.interface.h"

namespace saber {
namespace gsi {

// -----------------------------------------------------------------------------
UnstructuredInterpolation::UnstructuredInterpolation(const eckit::Configuration & config,
                                                     const atlas::FunctionSpace & fspace1,
                                                     const atlas::FunctionSpace & fspace2,
                                                     const std::vector<std::string> & activeVars,
                                                     const eckit::mpi::Comm & comm)
  : in_fspace_(&fspace1), out_fspace_(&fspace2), activeVars_(activeVars)
{
  saber_unstrc_create_f90(keyUnstructuredInterpolator_, &comm,
                          fspace1.lonlat().get(), fspace2.lonlat().get(), config);
}

// -----------------------------------------------------------------------------
int UnstructuredInterpolation::write(const eckit::Configuration & config) {
  saber_unstrc_write_f90(keyUnstructuredInterpolator_, config);
  return 0;
}

// -----------------------------------------------------------------------------
UnstructuredInterpolation::~UnstructuredInterpolation() {
  saber_unstrc_delete_f90(keyUnstructuredInterpolator_);
}
// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply(const atlas::Field & infield, atlas::Field & outfield) {
  saber_unstrc_apply_f90(keyUnstructuredInterpolator_, infield.get(), outfield.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply_ad(const atlas::Field & field_grid2,
                                         atlas::Field & field_grid1) {
  saber_unstrc_apply_ad_f90(keyUnstructuredInterpolator_, field_grid2.get(), field_grid1.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply(const atlas::FieldSet & innerFields,
                                      atlas::FieldSet & outerFields) {
  for (const auto & var : activeVars_) {
    ASSERT(innerFields.has(var));
    ASSERT(outerFields.has(var));
    this->apply(innerFields[var], outerFields[var]);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply_ad(const atlas::FieldSet & outerFields,
                                         atlas::FieldSet & innerFields) {
  for (const auto & var : activeVars_) {
    ASSERT(outerFields.has(var));
    ASSERT(innerFields.has(var));
    this->apply_ad(outerFields[var], innerFields[var]);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::print(std::ostream & os) const {
  os << " UnstructuredInterpolation: print not implemented yet.";
}
// -----------------------------------------------------------------------------
}  // namespace gsi
}  // namespace saber
