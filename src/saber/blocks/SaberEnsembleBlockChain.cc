/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberEnsembleBlockChain.h"

#include "oops/util/RandomField.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiply(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiply starting" << std::endl;

  if (strategy_ == "separated") {
    // Outer blocks adjoint multiplication
    if (outerBlockChain_) {
      outerBlockChain_->applyOuterBlocksAD(fset4d);
    }

    // Central block: ensemble covariance
    // Initialization
    const oops::FieldSet4D fset4dInit = oops::copyFieldSet4D(fset4d);
    fset4d.zero();

    for (const auto & scaleData : scaleDataVec_) {
      // Copy initial FieldSet4D
      oops::FieldSet4D fset4dScaleInit = oops::copyFieldSet4D(fset4dInit);

      if (scaleData.externalInterpolation()) {
        // Interpolate scale contribution to reduced resolution
        scaleData.interpolator()->applyOuterBlocksAD(fset4dScaleInit);
      }

      // Create scale FieldSet4D
      oops::FieldSet4D fset4dScale = oops::copyFieldSet4D(fset4dScaleInit);
      fset4dScale.zero();

      for (size_t ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
        // Copy initial FieldSet4D
        oops::FieldSet4D fset4dTmp = oops::copyFieldSet4D(fset4dScaleInit);

        if (scaleData.localization()) {
          // With localization

          // Get ensemble member
          oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
          if (scaleData.internalInterpolation()) {
            // Interpolate ensemble member to full resolution
            fset4dMem.deepCopy(*scaleData.ensemble(), ie);
            scaleData.interpolator()->applyOuterBlocks(fset4dMem);
          } else {
            // Keep ensemble member at file resolution
            for (size_t it = 0; it < fset4dMem.size(); ++it) {
              fset4dMem[it].shallowCopy((*scaleData.ensemble())(it, ie));
            }
          }

          // First schur product
          fset4dTmp *= fset4dMem;

          // Apply localization
          scaleData.localization()->multiply(fset4dTmp);

          // Second schur product
          fset4dTmp *= fset4dMem;

          // Add up member contribution
          fset4dScale += fset4dTmp;
        } else {
          // No localization

          // Compute weight
          const double wgt = fset4dScaleInit.dot_product_with(*scaleData.ensemble(), ie, vars_);

          // Copy ensemble member
          fset4dTmp.deepCopy(*scaleData.ensemble(), ie);

          // Apply weight
          fset4dTmp *= wgt;

          // Add up member contribution
          fset4dScale += fset4dTmp;
        }
        // TODO(Algo): Add communication here when the code starts supporting
        // ensemble members distributed across MPI tasks.
      }

      if (scaleData.externalInterpolation()) {
        // Interpolate scale contribution to full resolution
        scaleData.interpolator()->applyOuterBlocks(fset4dScale);
      }

      // Add up scale contribution
      fset4d += fset4dScale;
    }

    // Outer blocks forward multiplication
    if (outerBlockChain_) {
      outerBlockChain_->applyOuterBlocks(fset4d);
    }
  } else if (strategy_ == "crossed") {
    // Create control vector
    atlas::Field cv = atlas::Field("genericCtlVec",
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::make_shape(ctlVecSize()));

    // Adjoint square-root multiply (fset4D -> cv, index at 0)
    multiplySqrtAD(fset4d, cv, 0);

    // Square-root multiply (cv -> fset4D, index at 0)
    multiplySqrt(cv, fset4d, 0);
  }

  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::randomize(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "saber::SaberEnsembleBlockChain::randomize starting" << std::endl;

  if (strategy_ == "separated") {
    // Central block: randomization with ensemble covariance
    const auto & scaleData = scaleDataVec_[0];
    fset4d.deepCopy(*scaleData.ensemble(), 0);
    fset4d.zero();
    if (scaleData.internalInterpolation() || scaleData.externalInterpolation()) {
      // Prepare output FieldSet4D at full resolution
      scaleData.interpolator()->applyOuterBlocks(fset4d);
    }
    std::unique_ptr<util::NormalDistribution<double>> normalDist;

    for (const auto & scaleData : scaleDataVec_) {
      // Create scale FieldSet4D
      oops::FieldSet4D fset4dScale(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

      // Initialize scale contribution
      if (scaleData.externalInterpolation()) {
        // Process scale contribution at reduced resolution
        fset4dScale.deepCopy(*scaleData.ensemble(), 0);
      } else {
        // Process scale contribution at full resolution
        fset4dScale.deepCopy(fset4d);
      }
      fset4dScale.zero();

      for (unsigned int ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
        // Create empty FieldSet4D
        oops::FieldSet4D fset4dTmp(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

        if (scaleData.localization()) {
          // With localization

          // Randomize localization
          scaleData.localization()->randomize(fset4dTmp);

          // Get ensemble member
          oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
          if (scaleData.internalInterpolation()) {
            // Interpolate ensemble member to full resolution
            fset4dMem.deepCopy(*scaleData.ensemble(), ie);
            scaleData.interpolator()->applyOuterBlocks(fset4dMem);
          } else {
            // Keep ensemble member at file resolution
            for (size_t it = 0; it < fset4dMem.size(); ++it) {
              fset4dMem[it].shallowCopy((*scaleData.ensemble())(it, ie));
            }
          }

          // Schur product
          fset4dTmp *= fset4dMem;
        } else {
          // No localization
          if (!normalDist) {
            normalDist = std::make_unique<util::NormalDistribution<double>>(
              scaleData.ensemble()->ens_size(), 0.0, 1.0, seed_);
          }

          // Copy ensemble member
          fset4dTmp.deepCopy(*scaleData.ensemble(), ie);

          // Apply weight
          fset4dTmp *= (*normalDist)[ie];
        }

        // Add up member contribution
        fset4dScale += fset4dTmp;
      }

    if (scaleData.externalInterpolation()) {
        // Interpolate scale contribution to full resolution
        scaleData.interpolator()->applyOuterBlocks(fset4dScale);
      }

      // Add up scale contribution
      fset4d += fset4dScale;
    }

    // Outer blocks forward multiplication
    if (outerBlockChain_) {
       outerBlockChain_->applyOuterBlocks(fset4d);
    }
  } else if (strategy_ == "crossed") {
    // Create control vector
    atlas::Field cv("genericCtlVec", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(ctlVecSize()));

    // Generate random control vector
    randomCtlVec(cv, 0);

    // Square-root multiply
    multiplySqrt(cv, fset4d, 0);
  }

  oops::Log::trace() << "saber::SaberEnsembleBlockChain::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::randomCtlVec(atlas::Field & cv,
                                           const size_t & offset) const {
  oops::Log::trace() << "saber::SaberEnsembleBlockChain::randomCtlVec starting" << std::endl;

  // Initialization
  size_t index = offset;

  for (const auto & scaleData : scaleDataVec_) {
    if (strategy_ == "crossed") {
      // Restart index (save control vector)
      index = offset;
    }

    // Central block: ensemble covariance square-root
    for (unsigned int ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
      if (scaleData.localization()) {
        // With localization
        scaleData.localization()->randomCtlVec(cv, index);
        index += scaleData.localization()->ctlVecSize();
      } else {
        // No localization
        double zz;
        if (comm_.rank() == 0) {
          util::NormalDistributionField dist(1, 0.0, 1.0);
          zz = dist[0];
        }
        comm_.broadcast(zz, 0);
        auto cvView = atlas::array::make_view<double, 1>(cv);
        cvView(index) = zz;
        ++index;
      }
    }
  }

  oops::Log::trace() << "saber::SaberEnsembleBlockChain::randomCtlVec done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiplySqrt(const atlas::Field & cv,
                                           oops::FieldSet4D & fset4d,
                                           const size_t & offset) const {
  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiplySqrt starting" << std::endl;

  // Initialization
  const auto & scaleData = scaleDataVec_[0];
  fset4d.deepCopy(*scaleData.ensemble(), 0);
  fset4d.zero();
  if (scaleData.internalInterpolation() || scaleData.externalInterpolation()) {
    // Prepare output FieldSet4D at full resolution
    scaleData.interpolator()->applyOuterBlocks(fset4d);
  }
  size_t index = offset;

  for (const auto & scaleData : scaleDataVec_) {
    if (strategy_ == "crossed") {
      // Restart index (save control vector)
      index = offset;
    }

    // Create scale FieldSet4D
    oops::FieldSet4D fset4dScale(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

    // Initialize scale contribution
    if (scaleData.externalInterpolation()) {
      // Process scale contribution at reduced resolution
      fset4dScale.deepCopy(*scaleData.ensemble(), 0);
    } else {
      // Process scale contribution at full resolution
      fset4dScale.deepCopy(fset4d);
    }
    fset4dScale.zero();

    // Central block: ensemble covariance square-root
    for (unsigned int ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
      // Create empty FieldSet4D
      oops::FieldSet4D fset4dTmp(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());

      if (scaleData.localization()) {
        // With localization
        scaleData.localization()->multiplySqrt(cv, fset4dTmp, index);
        index += scaleData.localization()->ctlVecSize();

        // Get ensemble member
        oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
        if (scaleData.internalInterpolation()) {
          // Interpolate ensemble member to full resolution
          fset4dMem.deepCopy(*scaleData.ensemble(), ie);
          scaleData.interpolator()->applyOuterBlocks(fset4dMem);
        } else {
          // Keep ensemble member at file resolution
          for (size_t it = 0; it < fset4dMem.size(); ++it) {
            fset4dMem[it].shallowCopy((*scaleData.ensemble())(it, ie));
          }
        }

        // Schur product
        fset4dTmp *= fset4dMem;
      } else {
        // No localization
        const auto cvView = atlas::array::make_view<double, 1>(cv);

        // Copy ensemble member
        fset4dTmp.deepCopy(*scaleData.ensemble(), ie);

        // Apply weight
        fset4dTmp *= cvView(index);
        ++index;
      }

      // Add up member contribution
      fset4dScale += fset4dTmp;
    }

    if (scaleData.externalInterpolation()) {
      // Interpolate scale contribution to full resolution
      scaleData.interpolator()->applyOuterBlocks(fset4dScale);
    }

    // Add up scale contribution
    fset4d += fset4dScale;
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }

  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberEnsembleBlockChain::multiplySqrtAD(const oops::FieldSet4D & fset4d,
                                             atlas::Field & cv,
                                             const size_t & offset) const {
  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiplySqrtAD starting" << std::endl;

  // Copy input FieldSet
  oops::FieldSet4D fset4dInit = oops::copyFieldSet4D(fset4d);

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4dInit);
  }

  // Initialization
  size_t index = offset;

  // Get control vector view
  auto cvView = atlas::array::make_view<double, 1>(cv);

  // Initialize control vector
  cvView.assign(0.0);

  for (const auto & scaleData : scaleDataVec_) {
    if (strategy_ == "crossed") {
      // Restart index (save control vector)
      index = offset;
    }

    // Copy initial FieldSet4D
    oops::FieldSet4D fset4dScaleInit = oops::copyFieldSet4D(fset4dInit);

    if (scaleData.externalInterpolation()) {
      // Interpolate scale contribution to reduced resolution
      scaleData.interpolator()->applyOuterBlocksAD(fset4dScaleInit);
    }

    // Central block: ensemble covariance square-root adjoint
    if (scaleData.localization()) {
      // Create scale control vector
      atlas::Field cvScale("genericCtlVec", atlas::array::make_datatype<double>(),
        atlas::array::make_shape(scaleData.localization()->ctlVecSize()));

      // Get scale control vector view
      auto cvScaleView = atlas::array::make_view<double, 1>(cvScale);

      for (unsigned int ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
        // Copy initial fieldset
        oops::FieldSet4D fset4dTmp = oops::copyFieldSet4D(fset4dScaleInit);

        // Get ensemble member
        oops::FieldSet4D fset4dMem(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
        if (scaleData.internalInterpolation()) {
          // Interpolate ensemble member to full resolution
          fset4dMem.deepCopy(*scaleData.ensemble(), ie);
          scaleData.interpolator()->applyOuterBlocks(fset4dMem);
        } else {
          // Keep ensemble member at file resolution
          for (size_t it = 0; it < fset4dMem.size(); ++it) {
            fset4dMem[it].shallowCopy((*scaleData.ensemble())(it, ie));
          }
        }

        // Schur product
        fset4dTmp *= fset4dMem;

        // Apply localization square-root adjoint
        scaleData.localization()->multiplySqrtAD(fset4dTmp, cvScale, 0);
        for (size_t jcv = 0; jcv < scaleData.localization()->ctlVecSize(); ++jcv) {
          cvView(index+jcv) += cvScaleView(jcv);
        }
        index += scaleData.localization()->ctlVecSize();
      }
    } else {
      for (unsigned int ie = 0; ie < scaleData.ensemble()->ens_size(); ++ie) {
        // Compute weight
        cvView(index) = fset4dScaleInit.dot_product_with(*scaleData.ensemble(), ie, vars_);
        ++index;
      }
    }
  }

  oops::Log::trace() << "saber::SaberEnsembleBlockChain::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
