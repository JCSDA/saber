/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberEnsembleBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {


// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply starting" << std::endl;
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block: ensemble covariance
  // Initialization
  atlas::FieldSet fsetInit = util::copyFieldSet(fset);
  util::zeroFieldSet(fset);

  if (locBlockChain_) {
    // Apply localization
    for (size_t ie = 0; ie < ensemble_.size(); ++ie) {
      // Temporary copy for this ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(fsetInit);

      // First schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Apply localization
      locBlockChain_->multiply(fsetMem);

      // Second schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  } else {
    // No localization
    for (size_t ie = 0; ie < ensemble_.size(); ++ie) {
      // Compute weight
      const double wgt = util::dotProductFieldSets(fsetInit,
                                                   ensemble_[ie],
                                                   vars_.variables(),
                                                   comm_,
                                                   true);

      // Copy ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(ensemble_[ie]);

      // Apply weight
      util::multiplyFieldSet(fsetMem, wgt);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  }

  // Normalize result
  const double rk = 1.0/static_cast<double>(ensemble_.size()-1);
  util::multiplyFieldSet(fset, rk);

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }

  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::randomize(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block: randomization with ensemble covariance
  util::zeroFieldSet(fset);
  if (locBlockChain_) {
    // Apply localization
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Randomize SABER central block
      atlas::FieldSet fsetMem;
      locBlockChain_->randomize(fsetMem);

      // Schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  } else {
    // No localization
    util::NormalDistribution<double> normalDist(ensemble_.size(), 0.0, 1.0, seed_);
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Copy ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(ensemble_[ie]);

      // Apply weight
      util::multiplyFieldSet(fsetMem, normalDist[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  }

  // Normalize result
  const double rk = 1.0/sqrt(static_cast<double>(ensemble_.size()-1));
  util::multiplyFieldSet(fset, rk);

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

}  // namespace saber
