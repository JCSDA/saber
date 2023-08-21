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

void SaberParametricBlockChain::multiply(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block multiplication
  centralBlock_->multiply(fset);

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::randomize(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block randomization
  centralBlock_->randomize(fset);

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

}  // namespace saber
