/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberParametricBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiply(oops::FieldSet4D & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block multiplication
  for (size_t jtime = 0; jtime < fset.size(); ++jtime) {
    centralBlock_->multiply(fset[jtime].fieldSet());
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::randomize(oops::FieldSet4D & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block randomization
  for (size_t jtime = 0; jtime < fset.size(); ++jtime) {
    centralBlock_->randomize(fset[jtime].fieldSet());
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

}  // namespace saber
