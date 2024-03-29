/*
 * (C) Copyright 2024- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/GSIBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {

namespace gsi {

// -----------------------------------------------------------------------------

SaberGSIBlockChain::~SaberGSIBlockChain() {
  oops::Log::trace() << "SaberGSIBlockChain::~SaberGSIBlockChain starting" << std::endl;
  gsi_covariance_delete_f90(keySelf_);
  oops::Log::trace() << "SaberGSIBlockChain::~SaberGSIBlockChain done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberGSIBlockChain::multiply(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberGSIBlockChain::multiply starting" << std::endl;
  util::Timer timer("SaberGSIBlockChain", "multiply");

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4d);
  }

  // Central block multiplication
  std::vector<const atlas::field::FieldSetImpl*> fset4dptrs(fset4d.size());
  for (size_t itime = 0; itime < fset4d.size(); ++itime) {
    fset4dptrs[itime] = fset4d[itime].get();
  }
  gsi_covariance_multiply_f90(keySelf_, fset4dptrs.size(), fset4dptrs.data());

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
  oops::Log::trace() << "SaberGSIBlockChain::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------

void SaberGSIBlockChain::randomize(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberGSIBlockChain::randomize starting" << std::endl;
  util::Timer timer("SaberGSIBlockChain", "randomize");

  // Create central FieldSet4D
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    fset4d[jtime].init(centralFunctionSpace_, centralVars_);
  }

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on gsi grid
  for (const auto & sabField : fset4d[0]) {
      // Ensure that the field name is in the variables list
      if (!centralVars_.has(sabField.name())) {
        throw eckit::Exception("Field " + sabField.name() + " not found in the SaberGSIBlockChain" +
          " variables.", Here());
      }

      // Create the gsi grid field and add to Fieldset
      newFields.add(centralFunctionSpace_.createField<double>(atlas::option::name(sabField.name())
        | atlas::option::levels(sabField.shape(1))));
  }

  // Replace whatever fields are coming in with the gsi grid fields
  fset4d[0].fieldSet() = newFields;

  // Call implementation
  gsi_covariance_randomize_f90(keySelf_, fset4d[0].get());

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
  oops::Log::trace() << "SaberGSIBlockChain::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace gsi

}  // namespace saber
