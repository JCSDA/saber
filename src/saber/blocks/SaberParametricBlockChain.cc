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

void SaberParametricBlockChain::multiply(oops::FieldSet4D & fset4d) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4d);
  }

  // Central block multiplication
  if (crossTimeCov_) {
    // Duplicated cross-time covariances.
    // Use Mark Buehner's trick to save CPU when applying the same 3D cov/loc for all
    // 3D blocks of the 4D cov/loc matrix:
    // C_4D = ( C_3D C_3D C_3D ) = ( Id ) C_3D ( Id Id Id )
    //        ( C_3D C_3D C_3D )   ( Id )
    //        ( C_3D C_3D C_3D )   ( Id )
    // so if :
    // x_4D = ( x_1 )
    //        ( x_2 )
    //        ( x_3 )
    // then:
    // C_4D x_4D = (Id) C_3D (Id Id Id) (x_1) = (Id) C_3D (x_1+x_2+x_3) = (C_3D ( x_1 + x_2 + x_3 ))
    //             (Id)                 (x_2)   (Id)                      (C_3D ( x_1 + x_2 + x_3 ))
    //             (Id)                 (x_3)   (Id)                      (C_3D ( x_1 + x_2 + x_3 ))
    // Reference in section 3.4.2. of https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2325.
    // Local sum of x1, x2, ...
    for (size_t jtime = 1; jtime < fset4d.size(); ++jtime) {
      fset4d[0] += fset4d[jtime];
    }
    if (fset4d.commTime().rank() > 0) {
      oops::mpi::send(fset4d.commTime(), fset4d[0], 0, 0);
    } else {
      // On rank 0 receive other local sums and compute global sum of x1, x2, ...
      oops::FieldSet3D fset3d_tmp = oops::initFieldSet3D(fset4d[0]);
      for (size_t jj = 1; jj < fset4d.commTime().size(); ++jj) {
        oops::mpi::receive(fset4d.commTime(), fset3d_tmp, jj, 0);
        fset4d[0] += fset3d_tmp;
      }
      // Compute C * (x1+x2+...)
      centralBlock_->multiply(fset4d[0].fieldSet());
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(fset4d.commTime(), fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].fieldSet() = util::copyFieldSet(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances: apply central block to each of the
    // time slots.
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->multiply(fset4d[jtime].fieldSet());
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
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
