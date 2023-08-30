/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberEnsembleBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {


// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiply(oops::FieldSet4D & fset) const {
  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply starting" << std::endl;
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block: ensemble covariance
  // Initialization
  const oops::FieldSet4D fsetInit = oops::copyFieldSet4D(fset);
  fset.zero();

  if (locBlockChain_) {
    // Apply localization
    for (size_t ie = 0; ie < ensemble_.size(); ++ie) {
      // Temporary copy for this ensemble member
      oops::FieldSet4D fsetMem = oops::copyFieldSet4D(fsetInit);

      // First schur product
      fsetMem *= ensemble_[ie];

      // Apply localization
      locBlockChain_->multiply(fsetMem);

      // Second schur product
      fsetMem *= ensemble_[ie];

      // Add up member contribution
      fset += fsetMem;
    }
  } else {
    // No localization
    for (size_t ie = 0; ie < ensemble_.size(); ++ie) {
      // Compute weight
      const double wgt = util::dotProductFieldSets(fsetInit[0].fieldSet(),
                                                   ensemble_[ie],
                                                   vars_.variables(),
                                                   comm_,
                                                   true);

      // Copy ensemble member
      // TODO(AS): revisit ensemble_ to be a 4D ensemble, remove the hack below.
      oops::FieldSet4D fsetMem(fset.times(), fset.commTime(), fset[0].commGeom());
      for (size_t jtime = 0; jtime < fsetMem.size(); ++jtime) {
        fsetMem[jtime].fieldSet() = util::copyFieldSet(ensemble_[ie]);
      }

      // Apply weight
      fsetMem *= wgt;

      // Add up member contribution
      fset += fsetMem;
    }
  }

  // Normalize result
  const double rk = 1.0/static_cast<double>(ensemble_.size()-1);
  fset *= rk;

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }

  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::randomize(oops::FieldSet4D & fset) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset);
  }

  // Central block: randomization with ensemble covariance
  fset.zero();
  if (locBlockChain_) {
    // Apply localization
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Randomize SABER central block
      oops::FieldSet4D fsetMem(fset.times(), fset.commTime(), fset[0].commGeom());
      locBlockChain_->randomize(fsetMem);

      // Schur product
      fsetMem *= ensemble_[ie];

      // Add up member contribution
      fset += fsetMem;
    }
  } else {
    // No localization
    util::NormalDistribution<double> normalDist(ensemble_.size(), 0.0, 1.0, seed_);
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Copy ensemble member
      oops::FieldSet4D fsetMem(fset.times(), fset.commTime(), fset[0].commGeom());
      for (size_t jtime = 0; jtime < fsetMem.size(); ++jtime) {
        fsetMem[jtime].fieldSet() = util::copyFieldSet(ensemble_[ie]);
      }

      // Apply weight
      fsetMem *= normalDist[ie];

      // Add up member contribution
      fset += fsetMem;
    }
  }

  // Normalize result
  const double rk = 1.0/sqrt(static_cast<double>(ensemble_.size()-1));
  fset *= rk;

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset);
  }
}

}  // namespace saber
