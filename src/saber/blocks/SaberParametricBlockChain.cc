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

void SaberParametricBlockChain::filter(oops::FieldSet4D & fset4d) const {
  // Outer blocks for adjoint multiplication or left inverse (acting as filter)
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksFilter(fset4d);
  }

  // No cross-time covariances: apply central block to each of the
  // time slots.
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    centralBlock_->multiply(fset4d[jtime].fieldSet());
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

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
    if (timeComm_.rank() > 0) {
      oops::mpi::send(timeComm_, fset4d[0], 0, 0);
    } else {
      // On rank 0 receive other local sums and compute global sum of x1, x2, ...
      oops::FieldSet3D fset3d_tmp = oops::initFieldSet3D(fset4d[0]);
      for (size_t jj = 1; jj < timeComm_.size(); ++jj) {
        oops::mpi::receive(timeComm_, fset3d_tmp, jj, 0);
        fset4d[0] += fset3d_tmp;
      }
      // Compute C * (x1+x2+...)
      centralBlock_->multiply(fset4d[0].fieldSet());
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
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

void SaberParametricBlockChain::randomize(oops::FieldSet4D & fset4d) const {
  // Create central FieldSet4D
  atlas::FieldSet fset = util::createFieldSet(centralFunctionSpace_,
                                              centralVars_);
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    fset4d[jtime].fieldSet() = util::copyFieldSet(fset);
  }

  // Central block randomization
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      centralBlock_->randomize(fset4d[0].fieldSet());
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].fieldSet() = util::copyFieldSet(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->randomize(fset4d[jtime].fieldSet());
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

size_t SaberParametricBlockChain::ctlVecSize() const {
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      // Central block square-root for rank 0
      return centralBlock_->ctlVecSize();
    } else {
      // No control vector
      return 0;
    }
  } else {
    // No cross-time covariances
    return centralBlock_->ctlVecSize()*size4D_;
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiplySqrt(const atlas::Field & cv,
                                             oops::FieldSet4D & fset4d,
                                             const size_t & offset) const {
  // Create central FieldSet4D
  atlas::FieldSet fset = util::createFieldSet(centralFunctionSpace_,
                                              centralVars_);
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    fset4d[jtime].fieldSet() = util::copyFieldSet(fset);
  }

  // Central block square-root
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      // Central block square-root for rank 0
      centralBlock_->multiplySqrt(cv, fset4d[0].fieldSet(), offset);
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].fieldSet() = util::copyFieldSet(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances
    size_t index = offset;
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->multiplySqrt(cv, fset4d[jtime].fieldSet(), index);
      index += centralBlock_->ctlVecSize();
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiplySqrtAD(const oops::FieldSet4D & fset4d,
                                               atlas::Field & cv,
                                               const size_t & offset) const {
  // Initialization
  oops::FieldSet4D fset4dCopy = oops::copyFieldSet4D(fset4d);

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4dCopy);
  }

  // Central block square-root adjoint
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    for (size_t jtime = 1; jtime < fset4dCopy.size(); ++jtime) {
      fset4dCopy[0] += fset4dCopy[jtime];
    }
    if (timeComm_.rank() > 0) {
      oops::mpi::send(timeComm_, fset4dCopy[0], 0, 0);
    } else {
      // On rank 0 receive other local sums and compute global sum of x1, x2, ...
      oops::FieldSet3D fset3d_tmp = oops::initFieldSet3D(fset4dCopy[0]);
      for (size_t jj = 1; jj < timeComm_.size(); ++jj) {
        oops::mpi::receive(timeComm_, fset3d_tmp, jj, 0);
        fset4dCopy[0] += fset3d_tmp;
      }

      // Central block square-root adjoint for rank 0
      centralBlock_->multiplySqrtAD(fset4dCopy[0].fieldSet(), cv, offset);
    }
  } else {
    // No cross-time covariances
    size_t index = offset;
    for (size_t jtime = 0; jtime < fset4dCopy.size(); ++jtime) {
      centralBlock_->multiplySqrtAD(fset4dCopy[jtime].fieldSet(), cv, index);
      index += centralBlock_->ctlVecSize();
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
