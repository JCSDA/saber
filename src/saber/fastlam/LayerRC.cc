/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/LayerRC.h"

#include <algorithm>
#include <utility>

#include "atlas/array.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

static LayerMaker<LayerRC> makerRC_("rows-columns");

// -----------------------------------------------------------------------------

void LayerRC::setupParallelization() {
  oops::Log::trace() << classname() << "::setupParallelization starting" << std::endl;

  // Get index fields
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);

  // Sizes per task
  nxPerTask_.resize(comm_.size());
  nyPerTask_.resize(comm_.size());
  std::fill(nxPerTask_.begin(), nxPerTask_.end(), 0);
  std::fill(nyPerTask_.begin(), nyPerTask_.end(), 0);
  size_t index = 0;
  for (size_t jx = 0; jx < nx_; ++jx) {
    ++nxPerTask_[index];
    ++index;
    if (index == comm_.size()) index = 0;
  }
  index = 0;
  for (size_t jy = 0; jy < ny_; ++jy) {
    ++nyPerTask_[index];
    ++index;
    if (index == comm_.size()) index = 0;
  }

  // Start/end indices per task
  nxStart_.push_back(0);
  nyStart_.push_back(0);
  nxEnd_.push_back(nxPerTask_[0]-1);
  nyEnd_.push_back(nyPerTask_[0]-1);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    nxStart_.push_back(nxStart_[jt]+nxPerTask_[jt]);
    nyStart_.push_back(nyStart_[jt]+nyPerTask_[jt]);
    nxEnd_.push_back(nxStart_[jt+1]+nxPerTask_[jt+1]-1);
    nyEnd_.push_back(nyStart_[jt+1]+nyPerTask_[jt+1]-1);
  }

  // Rows <=> reduced grid

  // Define destination task
  xTask_.resize(rSize_);
  xOffset_.resize(rSize_);
  std::vector<int> xOffsetPerTask(comm_.size(), 0);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    bool found = false;
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      if (static_cast<size_t>(indexJView(jnode))-1 >= nyStart_[jt] &&
        static_cast<size_t>(indexJView(jnode))-1 <= nyEnd_[jt]) {
        xTask_[jnode] = jt;
        xOffset_[jnode] = xOffsetPerTask[jt];
        ++xOffsetPerTask[jt];
        ASSERT(!found);
        found = true;
      }
    }
    ASSERT(found);
  }

  // RecvCounts
  rRecvCounts_.resize(comm_.size());
  std::fill(rRecvCounts_.begin(), rRecvCounts_.end(), 0);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    ++rRecvCounts_[xTask_[jnode]];
  }

  // Buffer size
  rRecvSize_ = 0;
  for (const auto & n : rRecvCounts_) rRecvSize_ += n;

  // RecvDispls
  rRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    rRecvDispls_.push_back(rRecvDispls_[jt]+rRecvCounts_[jt]);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> rRecvCountsBuffer(comm_.size());
  comm_.allGatherv(rRecvCounts_.begin(), rRecvCounts_.end(), rRecvCountsBuffer);
  std::vector<int> rRecvCountsGlb_ = std::move(rRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xSendCounts_.push_back(rRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  xSendSize_ = 0;
  for (const auto & n : xSendCounts_) xSendSize_ += n;

  // SendDispls
  xSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xSendDispls_.push_back(xSendDispls_[jt]+xSendCounts_[jt]);
  }

  // Communicate indices
  std::vector<int> redIndex_i(rSize_);
  std::vector<int> redIndex_j(rSize_);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    redIndex_i[jnode] = indexIView(jnode)-1;
    redIndex_j[jnode] = indexJView(jnode)-1;
  }
  xIndex_i_.resize(xSendSize_);
  xIndex_j_.resize(xSendSize_);
  comm_.allToAllv(redIndex_i.data(), rRecvCounts_.data(), rRecvDispls_.data(),
    xIndex_i_.data(), xSendCounts_.data(), xSendDispls_.data());
  comm_.allToAllv(redIndex_j.data(), rRecvCounts_.data(), rRecvDispls_.data(),
    xIndex_j_.data(), xSendCounts_.data(), xSendDispls_.data());

  // Columns <=> rows

  // Define destination task
  yTask_.resize(nyPerTask_[myrank_]);
  yOffset_.resize(nyPerTask_[myrank_]);
  for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
    yTask_[j].resize(nx_);
    yOffset_[j].resize(nx_);
  }
  std::vector<int> yOffsetPerTask(comm_.size(), 0);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      bool found = false;
      for (size_t jt = 0; jt < comm_.size(); ++jt) {
        if (i >= nxStart_[jt] && i <= nxEnd_[jt]) {
          yTask_[j][i] = jt;
          yOffset_[j][i] = yOffsetPerTask[jt];
          ++yOffsetPerTask[jt];
          ASSERT(!found);
          found = true;
        }
      }
      ASSERT(found);
    }
  }

  // RecvCounts
  xRecvCounts_.resize(comm_.size());
  std::fill(xRecvCounts_.begin(), xRecvCounts_.end(), 0);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      ++xRecvCounts_[yTask_[j][i]];
    }
  }

  // Buffer size
  xRecvSize_ = 0;
  for (const auto & n : xRecvCounts_) xRecvSize_ += n;

  // RecvDispls
  xRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xRecvDispls_.push_back(xRecvDispls_[jt]+xRecvCounts_[jt]);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> xRecvCountsBuffer(comm_.size());
  comm_.allGatherv(xRecvCounts_.begin(), xRecvCounts_.end(), xRecvCountsBuffer);
  std::vector<int> xRecvCountsGlb_ = std::move(xRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ySendCounts_.push_back(xRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  ySendSize_ = 0;
  for (const auto & n : ySendCounts_) ySendSize_ += n;

  // SendDispls
  ySendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    ySendDispls_.push_back(ySendDispls_[jt]+ySendCounts_[jt]);
  }

  // Communicate indices
  std::vector<int> xIndex_i;
  std::vector<int> xIndex_j;
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      xIndex_i.push_back(i);
      xIndex_j.push_back(j+nyStart_[myrank_]);
    }
  }
  yIndex_i_.resize(ySendSize_);
  yIndex_j_.resize(ySendSize_);
  comm_.allToAllv(xIndex_i.data(), xRecvCounts_.data(), xRecvDispls_.data(),
    yIndex_i_.data(), ySendCounts_.data(), ySendDispls_.data());
  comm_.allToAllv(xIndex_j.data(), xRecvCounts_.data(), xRecvDispls_.data(),
    yIndex_j_.data(), ySendCounts_.data(), ySendDispls_.data());

  if (!params_.skipTests.value()) {
    // Tests

    // Generate fields
    atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy")
      | atlas::option::levels(nz_));
    atlas::Field rowsField = atlas::Field("dummy", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));
    atlas::Field colsField = atlas::Field("dummy", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
    auto redView = atlas::array::make_view<double, 2>(redField);
    auto rowsView = atlas::array::make_view<double, 3>(rowsField);
    auto colsView = atlas::array::make_view<double, 3>(colsField);

    // Fill field with global horizontal index
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        redView(jnode, k) = static_cast<double>(redIndex_i[jnode]*ny_+redIndex_j[jnode]);
      }
    }

    // Test rows <=> reduced grid
    redToRows(redField, rowsField);
    rowsToRed(rowsField, redField);

    // Print result
    oops::Log::test() << "    FastLAM redToRows test";
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (rowsView(i, j, k) != static_cast<double>(i*ny_+(j+nyStart_[myrank_]))) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("redToRows test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        if (redView(jnode, k) != static_cast<double>(redIndex_i[jnode]*ny_+redIndex_j[jnode])) {
          oops::Log::test() << " failed" << std::endl;
          throw eckit::Exception("redToRows test failed for block FastLAM", Here());
        }
      }
    }
    oops::Log::test() << " passed" << std::endl;

    // Test columns <=> rows
    rowsToCols(rowsField, colsField);
    colsToRows(colsField, rowsField);

    // Print result
    oops::Log::test() << "    FastLAM rowsToCols test";
    for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (colsView(i, j, k) != static_cast<double>((i+nxStart_[myrank_])*ny_+j)) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("rowsToCols test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (rowsView(i, j, k) != static_cast<double>(i*ny_+(j+nyStart_[myrank_]))) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("rowsToCols test failed for block FastLAM", Here());
          }
        }
      }
    }
    oops::Log::test() << " passed" << std::endl;
  }

  oops::Log::trace() << classname() << "::setupParallelization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::extractConvolution(const size_t & nxHalf,
                                 const size_t & nyHalf,
                                 std::vector<double> & horConv,
                                 std::vector<double> & verConv) {
  oops::Log::trace() << classname() << "::extractConvolution starting" << std::endl;

  // Reduced grid indices
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);

  // One level only
  const size_t nzSave = nz_;
  nz_ = 1;

  // Horizontal convolution
  atlas::Field redFieldHor = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));
  atlas::Field colsFieldHor("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto redViewHor = atlas::array::make_view<double, 2>(redFieldHor);
  for (size_t i = 0; i < nxHalf; ++i) {
    for (size_t j = 0; j < nyHalf; ++j) {
      // Setup dirac point
      redViewHor.assign(0.0);
      for (size_t jnode = 0; jnode < rSize_; ++jnode) {
        if (indexIView(jnode)-1 == static_cast<int>(i)
          && indexJView(jnode)-1 == static_cast<int>(j)) {
          redViewHor(jnode, 0) = 1.0;
        }
      }

      // Apply convolution on reduced grid
      multiplyRedSqrtTrans(redFieldHor, colsFieldHor);
      multiplyRedSqrt(colsFieldHor, redFieldHor);

      // Gather horizontal convolution data
      redViewHor = atlas::array::make_view<double, 2>(redFieldHor);
      for (size_t jnode = 0; jnode < rSize_; ++jnode) {
        if (indexIView(jnode)-1 == static_cast<int>(i)
          && indexJView(jnode)-1 == static_cast<int>(j)) {
          horConv[4*(i*nyHalf+j)+0] = redViewHor(jnode, 0);
        }
        if (indexIView(jnode)-1 == static_cast<int>(i+1)
          && indexJView(jnode)-1 == static_cast<int>(j)) {
          horConv[4*(i*nyHalf+j)+1] = redViewHor(jnode, 0);
        }
        if (indexIView(jnode)-1 == static_cast<int>(i)
          && indexJView(jnode)-1 == static_cast<int>(j+1)) {
          horConv[4*(i*nyHalf+j)+2] = redViewHor(jnode, 0);
        }
        if (indexIView(jnode)-1 == static_cast<int>(i+1)
          && indexJView(jnode)-1 == static_cast<int>(j+1)) {
          horConv[4*(i*nyHalf+j)+3] = redViewHor(jnode, 0);
        }
      }
    }
  }
  comm_.allReduceInPlace(horConv.begin(), horConv.end(), eckit::mpi::sum());

  // Reset number of levels
  nz_ = nzSave;

  // One grid point only
  const size_t nxPerTaskSave = nxPerTask_[myrank_];
  const size_t nySave = ny_;
  nxPerTask_[myrank_] = 1;
  ny_ = 1;

  // Vertical convolution
  atlas::Field colsVerField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto colsVerView = atlas::array::make_view<double, 3>(colsVerField);
  for (size_t k = 0; k < nz_-1; ++k) {
    // Setup dirac point
    colsVerView.assign(0.0);
    colsVerView(0, 0, k) = 1.0;

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Apply vertical kernel
    vertConvolution(colsVerField);
    vertConvolution(colsVerField);

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Gather horizontal convolution data
    verConv[2*k+0] = colsVerView(0, 0, k);
    verConv[2*k+1] = colsVerView(0, 0, k+1);
  }

  // Reset number of grid points
  nxPerTask_[myrank_] = nxPerTaskSave;
  ny_ = nySave;

  oops::Log::trace() << classname() << "::extractConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::multiplySqrt(const atlas::Field & cv,
                           atlas::Field & modelField,
                           const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Create field on columns
  atlas::Field colsField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));

  // Control vector to columns
  const auto cvView = atlas::array::make_view<double, 1>(cv);
  auto colsView = atlas::array::make_view<double, 3>(colsField);
  size_t index = offset;
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        colsView(i, j, k) = cvView(index);
        ++index;
      }
    }
  }

  // Create field on reduced grid
  atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));

  // Square-root multiplication on reduced grid
  multiplyRedSqrt(colsField, redField);

  // Interpolation TL
  interpolationTL(redField, modelField);

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::multiplySqrtTrans(const atlas::Field & modelField,
                                atlas::Field & cv,
                                const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtTrans starting" << std::endl;

  // Create field on reduced grid
  atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));

  // Interpolation AD
  interpolationAD(modelField, redField);

  // Create field on columns
  atlas::Field colsField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));

  // Adjoint square-root multiplication on reduced grid
  multiplyRedSqrtTrans(redField, colsField);

  // Columns to control vector
  const auto colsView = atlas::array::make_view<double, 3>(colsField);
  auto cvView = atlas::array::make_view<double, 1>(cv);
  size_t index = offset;
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        cvView(index) = colsView(i, j, k);
        ++index;
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

void LayerRC::redToRows(const atlas::Field & redField,
                        atlas::Field & rowsField) const {
  oops::Log::trace() << classname() << "::redToRows starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> rRecvCounts3D(comm_.size());
  std::vector<int> rRecvDispls3D(comm_.size());
  std::vector<int> xSendCounts3D(comm_.size());
  std::vector<int> xSendDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rRecvCounts3D[jt] = rRecvCounts_[jt]*nz_;
    rRecvDispls3D[jt] = rRecvDispls_[jt]*nz_;
    xSendCounts3D[jt] = xSendCounts_[jt]*nz_;
    xSendDispls3D[jt] = xSendDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> rRecvVec(rRecvSize_*nz_);
  const auto redView = atlas::array::make_view<double, 2>(redField);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = rRecvDispls3D[xTask_[jnode]] + xOffset_[jnode]*nz_ + k;
      rRecvVec[jv] = redView(jnode, k);
    }
  }

  // Communication
  std::vector<double> xSendVec(xSendSize_*nz_);
  comm_.allToAllv(rRecvVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data(),
    xSendVec.data(), xSendCounts3D.data(), xSendDispls3D.data());

  // Deserialize
  auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  for (size_t jx = 0; jx < xSendSize_; ++jx) {
    const size_t i = xIndex_i_[jx];
    const size_t j = xIndex_j_[jx]-nyStart_[myrank_];
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = jx*nz_ + k;
      rowsView(i, j, k) = xSendVec[jv];
    }
  }

  oops::Log::trace() << classname() << "::redToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::rowsToRed(const atlas::Field & rowsField,
                        atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::rowsToRed starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> rRecvCounts3D(comm_.size());
  std::vector<int> rRecvDispls3D(comm_.size());
  std::vector<int> xSendCounts3D(comm_.size());
  std::vector<int> xSendDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rRecvCounts3D[jt] = rRecvCounts_[jt]*nz_;
    rRecvDispls3D[jt] = rRecvDispls_[jt]*nz_;
    xSendCounts3D[jt] = xSendCounts_[jt]*nz_;
    xSendDispls3D[jt] = xSendDispls_[jt]*nz_;
  }

  // Serialize
  const auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  std::vector<double> xSendVec(xSendSize_*nz_);
  for (size_t jx = 0; jx < xSendSize_; ++jx) {
    const size_t i = xIndex_i_[jx];
    const size_t j = xIndex_j_[jx]-nyStart_[myrank_];
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = jx*nz_ + k;
      xSendVec[jv] = rowsView(i, j, k);
    }
  }

  // Communication
  std::vector<double> rRecvVec(rRecvSize_*nz_);
  comm_.allToAllv(xSendVec.data(), xSendCounts3D.data(), xSendDispls3D.data(),
    rRecvVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data());

  // Deserialize
  auto redView = atlas::array::make_view<double, 2>(redField);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = rRecvDispls3D[xTask_[jnode]] + xOffset_[jnode]*nz_ + k;
      redView(jnode, k) = rRecvVec[jv];
    }
  }

  oops::Log::trace() << classname() << "::rowsToRed done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::rowsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsConvolution starting" << std::endl;

  // Copy field
  atlas::Field copyField = field.clone();
  const auto copyView = atlas::array::make_view<double, 3>(copyField);

  // Apply kernel
  auto view = atlas::array::make_view<double, 3>(field);
  view.assign(0.0);
  for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t jk = 0; jk < xKernelSize_; ++jk) {
        const size_t ii = i-jk+(xKernelSize_-1)/2;
        if (ii >= 0 && ii < nx_) {
          for (size_t k = 0; k < nz_; ++k) {
            view(i, j, k) += copyView(ii, j, k)*xKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::rowsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsNormalization starting" << std::endl;

  // Apply normalization
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
    for (size_t i = 0; i < xNormSize_; ++i) {
      for (size_t k = 0; k < nz_; ++k) {
        view(i, j, k) *= xNorm_[i];
        view(nx_-1-i, j, k) *= xNorm_[i];
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::rowsToCols(const atlas::Field & rowsField,
                         atlas::Field & colsField) const {
  oops::Log::trace() << classname() << "::rowsToCols starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> xRecvCounts3D(comm_.size());
  std::vector<int> xRecvDispls3D(comm_.size());
  std::vector<int> ySendCounts3D(comm_.size());
  std::vector<int> ySendDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xRecvCounts3D[jt] = xRecvCounts_[jt]*nz_;
    xRecvDispls3D[jt] = xRecvDispls_[jt]*nz_;
    ySendCounts3D[jt] = ySendCounts_[jt]*nz_;
    ySendDispls3D[jt] = ySendDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> xRecvVec(xRecvSize_*nz_);
  const auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        const size_t jv = xRecvDispls3D[yTask_[j][i]] + yOffset_[j][i]*nz_ + k;
        xRecvVec[jv] = rowsView(i, j, k);
      }
    }
  }

  // Communication
  std::vector<double> ySendVec(ySendSize_*nz_);
  comm_.allToAllv(xRecvVec.data(), xRecvCounts3D.data(), xRecvDispls3D.data(),
    ySendVec.data(), ySendCounts3D.data(), ySendDispls3D.data());

  // Deserialize
  auto colsView = atlas::array::make_view<double, 3>(colsField);
  for (size_t jy = 0; jy < ySendSize_; ++jy) {
    const size_t i = yIndex_i_[jy]-nxStart_[myrank_];
    const size_t j = yIndex_j_[jy];
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = jy*nz_ + k;
      colsView(i, j, k) = ySendVec[jv];
    }
  }

  oops::Log::trace() << classname() << "::rowsToCols done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::colsToRows(const atlas::Field & colsField,
                         atlas::Field & rowsField) const {
  oops::Log::trace() << classname() << "::colsToRows starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> xRecvCounts3D(comm_.size());
  std::vector<int> xRecvDispls3D(comm_.size());
  std::vector<int> ySendCounts3D(comm_.size());
  std::vector<int> ySendDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xRecvCounts3D[jt] = xRecvCounts_[jt]*nz_;
    xRecvDispls3D[jt] = xRecvDispls_[jt]*nz_;
    ySendCounts3D[jt] = ySendCounts_[jt]*nz_;
    ySendDispls3D[jt] = ySendDispls_[jt]*nz_;
  }

  // Serialize
  const auto colsView = atlas::array::make_view<double, 3>(colsField);
  std::vector<double> ySendVec(ySendSize_*nz_);
  for (size_t jy = 0; jy < ySendSize_; ++jy) {
    size_t i = yIndex_i_[jy]-nxStart_[myrank_];
    size_t j = yIndex_j_[jy];
    for (size_t k = 0; k < nz_; ++k) {
      const size_t jv = jy*nz_ + k;
      ySendVec[jv] = colsView(i, j, k);
    }
  }

  // Communication
  std::vector<double> xRecvVec(xRecvSize_*nz_);
  comm_.allToAllv(ySendVec.data(), ySendCounts3D.data(), ySendDispls3D.data(),
    xRecvVec.data(), xRecvCounts3D.data(), xRecvDispls3D.data());

  // Deserialize
  auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        const size_t jv = xRecvDispls3D[yTask_[j][i]] + yOffset_[j][i]*nz_ + k;
        rowsView(i, j, k) = xRecvVec[jv];
      }
    }
  }

  oops::Log::trace() << classname() << "::colsToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::colsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsConvolution starting" << std::endl;

  // Copy field
  atlas::Field copyField = field.clone();
  const auto copyView = atlas::array::make_view<double, 3>(copyField);

  // Apply kernel
  auto view = atlas::array::make_view<double, 3>(field);
  view.assign(0.0);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t jk = 0; jk < yKernelSize_; ++jk) {
        const size_t jj = j-jk+(yKernelSize_-1)/2;
        if (jj >= 0 && jj < ny_) {
          for (size_t k = 0; k < nz_; ++k) {
            view(i, j, k) += copyView(i, jj, k)*yKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::colsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::colsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsNormalization starting" << std::endl;

  // Apply normalization
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < yNormSize_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        view(i, j, k) *= yNorm_[j];
        view(i, ny_-1-j, k) *= yNorm_[j];
      }
    }
  }

  oops::Log::trace() << classname() << "::colsNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::vertConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertConvolution starting" << std::endl;

  // Copy field
  atlas::Field copyField = field.clone();
  const auto copyView = atlas::array::make_view<double, 3>(copyField);

  // Apply kernel
  auto view = atlas::array::make_view<double, 3>(field);
  view.assign(0.0);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        for (size_t jk = 0; jk < zKernelSize_; ++jk) {
          const size_t kk = k-jk+(zKernelSize_-1)/2;
          if (kk >= 0 && kk < nz_) {
            view(i, j, k) += copyView(i, j, kk)*zKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::vertConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::vertNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertNormalization starting" << std::endl;

  // Apply normalization
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < zNormSize_; ++k) {
        view(i, j, k) *= zNorm_[k];
        view(i, j, nz_-1-k) *= zNorm_[k];
      }
    }
  }

  oops::Log::trace() << classname() << "::vertNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::multiplyRedSqrt(const atlas::Field & colsField,
                              atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrt starting" << std::endl;

  // Create intermediate fields
  atlas::Field colsFieldTmp = colsField.clone();
  atlas::Field rowsField("dummy",
                         atlas::array::make_datatype<double>(),
                         atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));

  if (nz_ > 1) {
    // Apply vertical kernel
    vertConvolution(colsFieldTmp);

    // Apply vertical normalization
    vertNormalization(colsFieldTmp);
  }

  // Apply kernel on columns
  colsConvolution(colsFieldTmp);

  // Apply normalization on columns
  colsNormalization(colsFieldTmp);

  // Columns to rows
  colsToRows(colsFieldTmp, rowsField);

  // Apply kernel on rows
  rowsConvolution(rowsField);

  // Apply normalization on rows
  rowsNormalization(rowsField);

  // Rows to reduced grid
  rowsToRed(rowsField, redField);

  oops::Log::trace() << classname() << "::multiplyRedSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerRC::multiplyRedSqrtTrans(const atlas::Field & redField,
                                   atlas::Field & colsField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans starting" << std::endl;

  // Create intermediate fields
  atlas::Field rowsField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));

  // Reduced grid to rows
  redToRows(redField, rowsField);

  // Apply normalization on rows
  rowsNormalization(rowsField);

  // Apply kernel on rows
  rowsConvolution(rowsField);

  // Rows to columns
  rowsToCols(rowsField, colsField);

  // Apply normalization on columns
  colsNormalization(colsField);

  // Apply kernel on columns
  colsConvolution(colsField);

  if (nz_ > 1) {
    // Apply vertical normalization
    vertNormalization(colsField);

    // Apply vertical kernel
    vertConvolution(colsField);
  }

  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
