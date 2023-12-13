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

void SaberEnsembleBlockChain::multiply(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply starting" << std::endl;

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4d);
  }

  // Central block: ensemble covariance
  // Initialization
  const oops::FieldSet4D fset4dInit = oops::copyFieldSet4D(fset4d);
  fset4d.zero();
  for (size_t ie = 0; ie < ensemble_.ens_size(); ++ie) {
    // Copy initial FieldSet4D
    oops::FieldSet4D fset4dMem = oops::copyFieldSet4D(fset4dInit);

    if (locBlockChain_) {
      // With localization
      // First schur product
      for (size_t it = 0; it < fset4dMem.size(); ++it) {
        fset4dMem[it] *= ensemble_(it, ie);
      }
      // Apply localization
      locBlockChain_->multiply(fset4dMem);
      // Second schur product
      for (size_t it = 0; it < fset4dMem.size(); ++it) {
        fset4dMem[it] *= ensemble_(it, ie);
      }
      // Add up member contribution
      fset4d += fset4dMem;
    } else {
      // No localization
      // Compute weight
      const double wgt = fset4dInit.dot_product_with(ensemble_, ie, vars_);
      // Copy ensemble member
      fset4dMem.deepCopy(ensemble_, ie);
      // Apply weight
      fset4dMem *= wgt;
      // Add up member contribution
      fset4d += fset4dMem;
    }
    // TODO(Algo): Add communication here when the code starts supporting
    // ensemble members distributed across MPI tasks.
  }

  // Normalize result
  const double rk = 1.0/static_cast<double>(ensemble_.ens_size()-1);
  fset4d *= rk;

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }

  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::randomize(oops::FieldSet4D & fset4d) const {
  // Central block: randomization with ensemble covariance
  fset4d.zero();
  std::unique_ptr<util::NormalDistribution<double>> normalDist;

  for (unsigned int ie = 0; ie < ensemble_.ens_size(); ++ie) {
    // Create empty FieldSet4D
    oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

    if (locBlockChain_) {
      // With localization

      // Randomize localization
      locBlockChain_->randomize(fset4dMem);

      // Schur product
      for (size_t it = 0; it < fset4dMem.size(); ++it) {
        fset4dMem[it] *= ensemble_(it, ie);
      }
    } else {
      // No localization
      if (!normalDist) {
        normalDist.reset(new util::NormalDistribution<double>(ensemble_.ens_size(),
                                                              0.0, 1.0, seed_));
      }
      // Copy ensemble member
      fset4dMem.deepCopy(ensemble_, ie);

      // Apply weight
      fset4dMem *= (*normalDist)[ie];
    }

    // Add up member contribution
    fset4d += fset4dMem;
  }

  // Normalize result
  const double rk = 1.0/sqrt(static_cast<double>(ensemble_.ens_size()-1));
  fset4d *= rk;

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiplySqrt(const atlas::Field & cv,
                                           oops::FieldSet4D & fset4d,
                                           const size_t & offset) const {
  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiplySqrt starting"
                     << std::endl;

  // Initialization
  fset4d.zero();
  size_t index = offset;

  // Central block: ensemble covariance square-root
  for (unsigned int ie = 0; ie < ensemble_.ens_size(); ++ie) {
    // Create empty FieldSet4D
    oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

    if (locBlockChain_) {
      // With localization
      locBlockChain_->multiplySqrt(cv, fset4dMem, index);
      index += locBlockChain_->ctlVecSize();

      // Schur product
      for (size_t it = 0; it < fset4dMem.size(); ++it) {
        fset4dMem[it] *= ensemble_(it, ie);
      }
    } else {
      // No localization
      const auto cvView = atlas::array::make_view<double, 1>(cv);
      fset4dMem.deepCopy(ensemble_, ie);

      // Apply weight
      fset4dMem *= cvView(index);
      ++index;
    }

    // Add up member contribution
    fset4d += fset4dMem;
  }

  // Normalize result
  const double rk = 1.0/std::sqrt(static_cast<double>(ensemble_.ens_size()-1));
  fset4d *= rk;

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }

  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiplySqrt done"
                     << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiplySqrtAD(const oops::FieldSet4D & fset4d,
                                             atlas::Field & cv,
                                             const size_t & offset) const {
  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiplySqrtAD starting"
                     << std::endl;

  // Copy input FieldSet
  oops::FieldSet4D fset4dInit = oops::copyFieldSet4D(fset4d);

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4dInit);
  }

  // Normalize initial fieldset
  const double rk = 1.0/std::sqrt(static_cast<double>(ensemble_.ens_size()-1));
  fset4dInit *= rk;

  // Initialization
  size_t index = offset;

  // Central block: ensemble covariance square-root adjoint
  for (unsigned int ie = 0; ie < ensemble_.ens_size(); ++ie) {
    if (locBlockChain_) {
      // Apply localization

      // Copy initial fieldset
      oops::FieldSet4D fset4dMem = oops::copyFieldSet4D(fset4dInit);

      // First schur product
      for (size_t it = 0; it < fset4dMem.size(); ++it) {
        fset4dMem[it] *= ensemble_(it, ie);
      }

      // Apply localization square-root adjoint
      locBlockChain_->multiplySqrtAD(fset4dMem, cv, index);
      index += locBlockChain_->ctlVecSize();
    } else {
      // No localization
      auto cvView = atlas::array::make_view<double, 1>(cv);

      // Compute weight
      cvView(index) = fset4dInit.dot_product_with(ensemble_, ie, vars_);
      ++index;
    }
  }

  oops::Log::trace() << "saber::generic::SaberEnsembleBlockChain::multiplySqrtAD done"
                     << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
