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

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.interface.h"

namespace saber {
namespace gsi {

// -----------------------------------------------------------------------------
UnstructuredInterpolation::UnstructuredInterpolation(
  const eckit::mpi::Comm & comm,
  const eckit::Configuration & config,
  const atlas::FunctionSpace & innerFuncSpace,
  const atlas::FunctionSpace & outerFuncSpace,
  const oops::Variables & activeVars)
  : innerFuncSpace_(innerFuncSpace), outerFuncSpace_(outerFuncSpace),
    activeVars_(activeVars)
{
  saber_unstrc_create_f90(keyUnstructuredInterpolator_, &comm,
                          innerFuncSpace_.lonlat().get(), outerFuncSpace_.lonlat().get(), config);
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
void UnstructuredInterpolation::apply(const atlas::Field & innerField,
                                      atlas::Field & outerField) {
  saber_unstrc_apply_f90(keyUnstructuredInterpolator_, innerField.get(), outerField.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::applyAD(const atlas::Field & outerField,
                                        atlas::Field & innerField) {
  saber_unstrc_apply_ad_f90(keyUnstructuredInterpolator_, outerField.get(), innerField.get());
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::apply(atlas::FieldSet & fset) {
  for (size_t i = 0; i < activeVars_.size(); ++i) {
    const size_t levels = fset[activeVars_[i].name()].levels();
    atlas::Field outerField = outerFuncSpace_.createField<double>(
      atlas::option::name(activeVars_[i].name()) | atlas::option::levels(levels));
    this->apply(fset[activeVars_[i].name()], outerField);
    util::removeFieldsFromFieldSet(fset, {activeVars_[i].name()});
    fset.add(outerField);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::applyAD(atlas::FieldSet & fset) {
  for (size_t i = 0; i < activeVars_.size(); ++i) {
    const size_t levels = fset[activeVars_[i].name()].levels();
    atlas::Field innerField = innerFuncSpace_.createField<double>(
      atlas::option::name(activeVars_[i].name()) | atlas::option::levels(levels));
    this->applyAD(fset[activeVars_[i].name()], innerField);
    util::removeFieldsFromFieldSet(fset, {activeVars_[i].name()});
    fset.add(innerField);
  }
}

// -----------------------------------------------------------------------------
void UnstructuredInterpolation::print(std::ostream & os) const {
  os << " UnstructuredInterpolation: print not implemented yet.";
}
// -----------------------------------------------------------------------------
}  // namespace gsi
}  // namespace saber
