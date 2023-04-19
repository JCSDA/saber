/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberBlockChain.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberBlockChain::SaberBlockChain(const oops::Variables & incVars, const atlas::FieldSet & fset) :
  centralBlock_(), outerBlocks_(), centralFieldSet_(), incVars_(incVars), scalarWeightSqrt_(1.0),
  fileWeightSqrt_() {
  // Initialize central FieldSet
  centralFieldSet_ = util::copyFieldSet(fset);
  util::zeroFieldSet(centralFieldSet_);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::centralBlockInit(SaberCentralBlockBase * centralBlockPtr) {
  centralBlock_.reset(centralBlockPtr);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::setWeight(const double & scalarWeight) {
  ASSERT(scalarWeight > 0.0);
  scalarWeightSqrt_ = std::sqrt(scalarWeight);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::setWeight(const atlas::FieldSet & fileWeight) {
  fileWeightSqrt_ = util::copyFieldSet(fileWeight);
  util::sqrtFieldSet(fileWeightSqrt_);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::applyOuterBlocks(atlas::FieldSet & fset) const {
  // Outer blocks forward multiplication
  for (ircst_ it = outerBlocks_.rbegin(); it != outerBlocks_.rend(); ++it) {
    it->multiply(fset);
  }

  // Weight square-root multiplication
  if (fileWeightSqrt_.empty()) {
    // Scalar weight
    util::multiplyFieldSet(fset, scalarWeightSqrt_);
  } else {
    // File-based weight
    util::multiplyFieldSets(fset, fileWeightSqrt_);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::applyOuterBlocksAD(atlas::FieldSet & fset) const {
  // Weight square-root multiplication
  if (fileWeightSqrt_.empty()) {
    // Scalar weight
    util::multiplyFieldSet(fset, scalarWeightSqrt_);
  } else {
    // File-based weight
    util::multiplyFieldSets(fset, fileWeightSqrt_);
  }

  // Outer blocks adjoint multiplication
  for (icst_ it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
    it->multiplyAD(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::randomize(atlas::FieldSet & fset) const {
  // Initialize dx.fieldSet() with all the required fields, set to zero
  fset = util::copyFieldSet(centralFieldSet_);

  // Central block randomization
  centralBlock_->randomize(fset);

  // Outer blocks forward multiplication
  applyOuterBlocks(fset);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::multiply(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  applyOuterBlocksAD(fset);

  // Central block multiplication
  centralBlock_->multiply(fset);

  // Outer blocks forward multiplication
  applyOuterBlocks(fset);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::leftInverseMultiply(atlas::FieldSet & fset) const {
  // Outer blocks calibration inverse multiplication
  for (icst_ it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
    it->leftInverseMultiply(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::leftInverseMultiplyExceptLast(atlas::FieldSet & fset) const {
  // Outer blocks calibration inverse multiplication
  for (icst_ it = outerBlocks_.begin(); it != std::prev(outerBlocks_.end()); ++it) {
    it->leftInverseMultiply(fset);
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
