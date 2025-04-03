/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/LayerSpec.h"

#include <algorithm>
#include <utility>

#include "atlas/array.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

static LayerMaker<LayerSpec> makerSpec_("spectral");

// -----------------------------------------------------------------------------

void LayerSpec::setupParallelization() {
  oops::Log::trace() << classname() << "::setupParallelization starting" << std::endl;

  // Get index fields
  atlas::Field indexXField = fset_["indexX"];
  atlas::Field indexYField = fset_["indexY"];
  auto indexXView = atlas::array::make_view<int, 1>(indexXField);
  auto indexYView = atlas::array::make_view<int, 1>(indexYField);

  // Extended sizes
  nxExt_ = nx_+xKernelSize_;
  nyExt_ = ny_+yKernelSize_;

  // Spectral sizes
  nk_ = nxExt_/2+1;
  nl_ = nyExt_/2+1;

  // Sizes per task
  nxPerTask_.resize(comm_.size());
  nyPerTask_.resize(comm_.size());
  std::fill(nxPerTask_.begin(), nxPerTask_.end(), 0);
  std::fill(nyPerTask_.begin(), nyPerTask_.end(), 0);
  size_t index = 0;
  for (size_t jx = 0; jx < nxExt_; ++jx) {
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
  xSendTask_.resize(rSize_);
  xSendOffset_.resize(rSize_);
  std::vector<int> xSendOffsetPerTask(comm_.size(), 0);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    bool found = false;
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      if (static_cast<size_t>(indexYView(jnode))-1 >= nyStart_[jt] &&
        static_cast<size_t>(indexYView(jnode))-1 <= nyEnd_[jt]) {
        xSendTask_[jnode] = jt;
        xSendOffset_[jnode] = xSendOffsetPerTask[jt];
        ++xSendOffsetPerTask[jt];
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
    ++rRecvCounts_[xSendTask_[jnode]];
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
  std::vector<int> redIndex_x(rSize_);
  std::vector<int> redIndex_y(rSize_);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    redIndex_x[jnode] = indexXView(jnode)-1;
    redIndex_y[jnode] = indexYView(jnode)-1;
  }
  xSendIndex_x_.resize(xSendSize_);
  xSendIndex_y_.resize(xSendSize_);
  comm_.allToAllv(redIndex_x.data(), rRecvCounts_.data(), rRecvDispls_.data(),
    xSendIndex_x_.data(), xSendCounts_.data(), xSendDispls_.data());
  comm_.allToAllv(redIndex_y.data(), rRecvCounts_.data(), rRecvDispls_.data(),
    xSendIndex_y_.data(), xSendCounts_.data(), xSendDispls_.data());

  // Columns <=> rows

  // Define destination task
  ySendTask_.resize(nxExt_*nyPerTask_[myrank_]);
  ySendOffset_.resize(nxExt_*nyPerTask_[myrank_]);
  std::vector<int> ySendOffsetPerTask(comm_.size(), 0);
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      bool found = false;
      for (size_t jt = 0; jt < comm_.size(); ++jt) {
        if (jx >= nxStart_[jt] && jx <= nxEnd_[jt]) {
          ySendTask_[jx*nyPerTask_[myrank_]+jy] = jt;
          ySendOffset_[jx*nyPerTask_[myrank_]+jy] = ySendOffsetPerTask[jt];
          ++ySendOffsetPerTask[jt];
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
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      ++xRecvCounts_[ySendTask_[jx*nyPerTask_[myrank_]+jy]];
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
  std::vector<int> xSendIndex_x(nxExt_*nyPerTask_[myrank_]);
  std::vector<int> xSendIndex_y(nxExt_*nyPerTask_[myrank_]);
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      xSendIndex_x[jx*nyPerTask_[myrank_]+jy] = jx;
      xSendIndex_y[jx*nyPerTask_[myrank_]+jy] = jy+nyStart_[myrank_];
    }
  }
  ySendIndex_x_.resize(ySendSize_);
  ySendIndex_y_.resize(ySendSize_);
  comm_.allToAllv(xSendIndex_x.data(), xRecvCounts_.data(), xRecvDispls_.data(),
    ySendIndex_x_.data(), ySendCounts_.data(), ySendDispls_.data());
  comm_.allToAllv(xSendIndex_y.data(), xRecvCounts_.data(), xRecvDispls_.data(),
    ySendIndex_y_.data(), ySendCounts_.data(), ySendDispls_.data());

  // Rows FFTW setup
  int xRank = 1;
  int xN[] = {static_cast<int>(nxExt_)};
  int xHowmany = nyPerTask_[myrank_]*nz_;
  int *xInembed = NULL;
  const int xIstride = 1;
  const int xIdist = static_cast<int>(nxExt_);
  int *xOnembed = NULL;
  const int xOstride = 1;
  const int xOdist = static_cast<int>(nk_);
  xBufR_ = fftw_alloc_real(nxExt_*nyPerTask_[myrank_]*nz_);
  xBufC_ = fftw_alloc_complex(nk_*nyPerTask_[myrank_]*nz_);
  xPlan_r2c_ = fftw_plan_many_dft_r2c(xRank, xN, xHowmany, xBufR_, xInembed, xIstride, xIdist,
    xBufC_, xOnembed, xOstride, xOdist, FFTW_PATIENT);
  xPlan_c2r_ = fftw_plan_many_dft_c2r(xRank, xN, xHowmany, xBufC_, xOnembed, xOstride, xOdist,
    xBufR_, xInembed, xIstride, xIdist, FFTW_PATIENT);

  // Rows normalization factor
  xNormFFT_ = 1.0/static_cast<double>(nxExt_);

  // Rows spectral standard deviation
  double *xBufR1d = fftw_alloc_real(nxExt_);
  fftw_complex *xBufC1d = fftw_alloc_complex(nk_);
  fftw_plan xPlan_r2c1d = fftw_plan_dft_r2c_1d(nxExt_, xBufR1d, xBufC1d, FFTW_PATIENT);
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    if (jx <= (xKernelSize_-1)/2) {
      xBufR1d[jx] = xKernel_[(xKernelSize_-1)/2+jx];
    } else if (jx > nxExt_-1-(xKernelSize_-1)/2) {
      xBufR1d[jx] = xKernel_[jx-(nxExt_-(xKernelSize_-1)/2)];
    } else {
      xBufR1d[jx] = 0.0;
    }
  }
  fftw_execute(xPlan_r2c1d);
  for (size_t jk = 0; jk < nk_; ++jk) {
    xSpecStdDev_.push_back(xBufC1d[jk][0]);
  }

  // Columns FFTW setup
  int yRank = 1;
  int yN[] = {static_cast<int>(nyExt_)};
  int yHowmany = nxPerTask_[myrank_]*nz_;
  int *yInembed = NULL;
  const int yIstride = 1;
  const int yIdist = static_cast<int>(nyExt_);
  int *yOnembed = NULL;
  const int yOstride = 1;
  const int yOdist = static_cast<int>(nl_);
  yBufR_ = fftw_alloc_real(nxPerTask_[myrank_]*nyExt_*nz_);
  yBufC_ = fftw_alloc_complex(nxPerTask_[myrank_]*nl_*nz_);
  yPlan_r2c_ = fftw_plan_many_dft_r2c(yRank, yN, yHowmany, yBufR_, yInembed, yIstride, yIdist,
    yBufC_, yOnembed, yOstride, yOdist, FFTW_PATIENT);
  yPlan_c2r_ = fftw_plan_many_dft_c2r(yRank, yN, yHowmany, yBufC_, yOnembed, yOstride, yOdist,
    yBufR_, yInembed, yIstride, yIdist, FFTW_PATIENT);

  // Columns normalization factor
  yNormFFT_ = 1.0/static_cast<double>(nyExt_);

  // Rows spectral standard deviation
  double *yBufR1d = fftw_alloc_real(nyExt_);
  fftw_complex *yBufC1d = fftw_alloc_complex(nl_);
  fftw_plan yPlan_r2c1d = fftw_plan_dft_r2c_1d(nyExt_, yBufR1d, yBufC1d, FFTW_PATIENT);
  for (size_t jy = 0; jy < nyExt_; ++jy) {
    if (jy <= (yKernelSize_-1)/2) {
      yBufR1d[jy] = yKernel_[(yKernelSize_-1)/2+jy];
    } else if (jy > nyExt_-1-(yKernelSize_-1)/2) {
      yBufR1d[jy] = yKernel_[jy-(nyExt_-(yKernelSize_-1)/2)];
    } else {
      yBufR1d[jy] = 0.0;
    }
  }
  fftw_execute(yPlan_r2c1d);
  for (size_t jl = 0; jl < nl_; ++jl) {
    ySpecStdDev_.push_back(yBufC1d[jl][0]);
  }

  if (!params_.skipTests.value()) {
    // Tests

    // Generate fields
    atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy")
      | atlas::option::levels(nz_));
    atlas::Field rowsField = atlas::Field("dummy", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nxExt_, nyPerTask_[myrank_], nz_));
    atlas::Field colsField = atlas::Field("dummy", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nxPerTask_[myrank_], nyExt_, nz_));
    auto redView = atlas::array::make_view<double, 2>(redField);
    auto rowsView = atlas::array::make_view<double, 3>(rowsField);
    auto colsView = atlas::array::make_view<double, 3>(colsField);

    // Fill field with global horizontal index
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        redView(jnode, jz) = static_cast<double>(redIndex_x[jnode]*ny_+redIndex_y[jnode]);
      }
    }

    // Test rows <=> reduced grid
    redToRows(redField, rowsField);
    rowsToRed(rowsField, redField);

    // Print result
    oops::Log::test() << "    FastLAM redToRows test";
    for (size_t jx = 0; jx < nxExt_; ++jx) {
      for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
        for (size_t jz = 0; jz < nz_; ++jz) {
          double testValue = (jx < nx_) ? static_cast<double>(jx*ny_+(jy+nyStart_[myrank_])) : 0.0;
          if (rowsView(jx, jy, jz) != testValue) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("redToRows test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        if (redView(jnode, jz) != static_cast<double>(redIndex_x[jnode]*ny_+redIndex_y[jnode])) {
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
    for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
      for (size_t jy = 0; jy < nyExt_; ++jy) {
        for (size_t jz = 0; jz < nz_; ++jz) {
          double testValue = (jx+nxStart_[myrank_] < nx_ && jy < ny_) ?
            static_cast<double>((jx+nxStart_[myrank_])*ny_+jy) : 0.0;
          if (colsView(jx, jy, jz) != testValue) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("rowsToCols test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t jx = 0; jx < nx_; ++jx) {
      for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
        for (size_t jz = 0; jz < nz_; ++jz) {
          double testValue = (jx < nx_) ? static_cast<double>(jx*ny_+(jy+nyStart_[myrank_])) : 0.0;
          if (rowsView(jx, jy, jz) != testValue) {
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

void LayerSpec::extractConvolution(const size_t & nxHalf,
                                   const size_t & nyHalf,
                                   std::vector<double> & horConv,
                                   std::vector<double> & verConv) {
  oops::Log::trace() << classname() << "::extractConvolution starting" << std::endl;

  // Reduced grid indices
  atlas::Field indexXField = fset_["indexX"];
  atlas::Field indexYField = fset_["indexY"];
  auto indexXView = atlas::array::make_view<int, 1>(indexXField);
  auto indexYView = atlas::array::make_view<int, 1>(indexYField);

  // One level only
  const size_t nzSave = nz_;
  nz_ = 1;

  // Horizontal convolution
  atlas::Field redFieldHor = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));
  atlas::Field colsFieldHor("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], nyExt_, nz_));
  auto redViewHor = atlas::array::make_view<double, 2>(redFieldHor);
  for (size_t jx = 0; jx < nxHalf; ++jx) {
    for (size_t jy = 0; jy < nyHalf; ++jy) {
      // Setup dirac point
      redViewHor.assign(0.0);
      for (size_t jnode = 0; jnode < rSize_; ++jnode) {
        if (indexXView(jnode)-1 == static_cast<int>(jx)
          && indexYView(jnode)-1 == static_cast<int>(jy)) {
          redViewHor(jnode, 0) = 1.0;
        }
      }

      // Apply convolution on reduced grid
      multiplyRedSqrtTrans(redFieldHor, colsFieldHor);
      multiplyRedSqrt(colsFieldHor, redFieldHor);

      // Gather horizontal convolution data
      redViewHor = atlas::array::make_view<double, 2>(redFieldHor);
      for (size_t jnode = 0; jnode < rSize_; ++jnode) {
        if (indexXView(jnode)-1 == static_cast<int>(jx)
          && indexYView(jnode)-1 == static_cast<int>(jy)) {
          horConv[4*(jx*nyHalf+jy)+0] = redViewHor(jnode, 0);
        }
        if (indexXView(jnode)-1 == static_cast<int>(jx+1)
          && indexYView(jnode)-1 == static_cast<int>(jy)) {
          horConv[4*(jx*nyHalf+jy)+1] = redViewHor(jnode, 0);
        }
        if (indexXView(jnode)-1 == static_cast<int>(jx)
          && indexYView(jnode)-1 == static_cast<int>(jy+1)) {
          horConv[4*(jx*nyHalf+jy)+2] = redViewHor(jnode, 0);
        }
        if (indexXView(jnode)-1 == static_cast<int>(jx+1)
          && indexYView(jnode)-1 == static_cast<int>(jy+1)) {
          horConv[4*(jx*nyHalf+jy)+3] = redViewHor(jnode, 0);
        }
      }
    }
  }
  comm_.allReduceInPlace(horConv.begin(), horConv.end(), eckit::mpi::sum());

  // Reset number of levels
  nz_ = nzSave;

  // One grid point only
  const size_t nxPerTaskSave = nxPerTask_[myrank_];
  const size_t nyExtSave = nyExt_;
  nxPerTask_[myrank_] = 1;
  nyExt_ = 1;

  // Vertical convolution
  atlas::Field colsVerField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], nyExt_, nz_));
  auto colsVerView = atlas::array::make_view<double, 3>(colsVerField);
  for (size_t jz = 0; jz < nz_-1; ++jz) {
    // Setup dirac point
    colsVerView.assign(0.0);
    colsVerView(0, 0, jz) = 1.0;

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Apply vertical kernel
    vertConvolution(colsVerField);
    vertConvolution(colsVerField);

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Gather horizontal convolution data
    verConv[2*jz+0] = colsVerView(0, 0, jz);
    verConv[2*jz+1] = colsVerView(0, 0, jz+1);
  }

  // Reset number of grid points
  nxPerTask_[myrank_] = nxPerTaskSave;
  nyExt_ = nyExtSave;

  oops::Log::trace() << classname() << "::extractConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::multiplySqrt(const atlas::Field & cv,
                             atlas::Field & modelField,
                             const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Create field on columns
  atlas::Field colsField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], nyExt_, nz_));

  // Control vector to columns
  const auto cvView = atlas::array::make_view<double, 1>(cv);
  auto colsView = atlas::array::make_view<double, 3>(colsField);
  size_t index = offset;
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jy = 0; jy < nyExt_; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        colsView(jx, jy, jz) = cvView(index);
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

void LayerSpec::multiplySqrtTrans(const atlas::Field & modelField,
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
    atlas::array::make_shape(nxPerTask_[myrank_], nyExt_, nz_));

  // Adjoint square-root multiplication on reduced grid
  multiplyRedSqrtTrans(redField, colsField);

  // Columns to control vector
  const auto colsView = atlas::array::make_view<double, 3>(colsField);
  auto cvView = atlas::array::make_view<double, 1>(cv);
  size_t index = offset;
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jy = 0; jy < nyExt_; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        cvView(index) = colsView(jx, jy, jz);
        ++index;
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

void LayerSpec::redToRows(const atlas::Field & redField,
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
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = rRecvDispls3D[xSendTask_[jnode]] + xSendOffset_[jnode]*nz_ + jz;
      rRecvVec[jv] = redView(jnode, jz);
    }
  }

  // Communication
  std::vector<double> xSendVec(xSendSize_*nz_);
  comm_.allToAllv(rRecvVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data(),
    xSendVec.data(), xSendCounts3D.data(), xSendDispls3D.data());

  // Deserialize
  auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  for (size_t js = 0; js < xSendSize_; ++js) {
    size_t jx = xSendIndex_x_[js];
    size_t jy = xSendIndex_y_[js]-nyStart_[myrank_];
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = js*nz_ + jz;
      rowsView(jx, jy, jz) = xSendVec[jv];
    }
  }

  // Extend
  for (size_t jx = nx_; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        rowsView(jx, jy, jz) = 0.0;
      }
    }
  }

  oops::Log::trace() << classname() << "::redToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::rowsToRed(const atlas::Field & rowsField,
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
  for (size_t js = 0; js < xSendSize_; ++js) {
    const size_t jx = xSendIndex_x_[js];
    const size_t jy = xSendIndex_y_[js]-nyStart_[myrank_];
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = js*nz_ + jz;
      xSendVec[jv] = rowsView(jx, jy, jz);
    }
  }

  // Communication
  std::vector<double> rRecvVec(rRecvSize_*nz_);
  comm_.allToAllv(xSendVec.data(), xSendCounts3D.data(), xSendDispls3D.data(),
    rRecvVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data());

  // Deserialize
  auto redView = atlas::array::make_view<double, 2>(redField);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = rRecvDispls3D[xSendTask_[jnode]] + xSendOffset_[jnode]*nz_ + jz;
      redView(jnode, jz) = rRecvVec[jv];
    }
  }

  oops::Log::trace() << classname() << "::rowsToRed done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::rowsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsConvolution starting" << std::endl;

  // Serialize
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jx = 0; jx < nxExt_; ++jx) {
        const size_t jv = jy*(nz_*nxExt_) + jz*nxExt_ + jx;
        xBufR_[jv] = view(jx, jy, jz);
      }
    }
  }

  // Compute direct transform
  fftw_execute(xPlan_r2c_);

  // Convolution in spectral space
  for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jk = 0; jk < nk_; ++jk) {
        const size_t jv = jy*(nz_*nk_) + jz*nk_ + jk;
        xBufC_[jv][0] *= xSpecStdDev_[jk];
        xBufC_[jv][1] *= xSpecStdDev_[jk];
      }
    }
  }

  // Compute inverse transform
  fftw_execute(xPlan_c2r_);

  // Deserialize and normalize
  for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jx = 0; jx < nxExt_; ++jx) {
        const size_t jv = jy*(nz_*nxExt_) + jz*nxExt_ + jx;
        view(jx, jy, jz) = xBufR_[jv]*xNormFFT_;
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::rowsToCols(const atlas::Field & rowsField,
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
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        const size_t jv = xRecvDispls3D[ySendTask_[jx*nyPerTask_[myrank_]+jy]]
          + ySendOffset_[jx*nyPerTask_[myrank_]+jy]*nz_ + jz;
        xRecvVec[jv] = rowsView(jx, jy, jz);
      }
    }
  }

  // Communication
  std::vector<double> ySendVec(ySendSize_*nz_);
  comm_.allToAllv(xRecvVec.data(), xRecvCounts3D.data(), xRecvDispls3D.data(),
    ySendVec.data(), ySendCounts3D.data(), ySendDispls3D.data());

  // Deserialize
  auto colsView = atlas::array::make_view<double, 3>(colsField);
  for (size_t js = 0; js < ySendSize_; ++js) {
    const size_t jx = ySendIndex_x_[js]-nxStart_[myrank_];
    const size_t jy = ySendIndex_y_[js];
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = js*nz_ + jz;
      colsView(jx, jy, jz) = ySendVec[jv];
    }
  }

  // Extend
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jy = ny_; jy < nyExt_; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        colsView(jx, jy, jz) = 0.0;
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsToCols done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::colsToRows(const atlas::Field & colsField,
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
  for (size_t js = 0; js < ySendSize_; ++js) {
    const size_t jx = ySendIndex_x_[js]-nxStart_[myrank_];
    const size_t jy = ySendIndex_y_[js];
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jv = js*nz_ + jz;
      ySendVec[jv] = colsView(jx, jy, jz);
    }
  }

  // Communication
  std::vector<double> xRecvVec(xRecvSize_*nz_);
  comm_.allToAllv(ySendVec.data(), ySendCounts3D.data(), ySendDispls3D.data(),
    xRecvVec.data(), xRecvCounts3D.data(), xRecvDispls3D.data());

  // Deserialize
  auto rowsView = atlas::array::make_view<double, 3>(rowsField);
  for (size_t jx = 0; jx < nxExt_; ++jx) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        const size_t jv = xRecvDispls3D[ySendTask_[jx*nyPerTask_[myrank_]+jy]]
          + ySendOffset_[jx*nyPerTask_[myrank_]+jy]*nz_ + jz;
        rowsView(jx, jy, jz) = xRecvVec[jv];
      }
    }
  }

  oops::Log::trace() << classname() << "::colsToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::colsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsConvolution starting" << std::endl;

  // Serialize
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jy = 0; jy < nyExt_; ++jy) {
        const size_t jv = jx*(nz_*nyExt_) + jz*nyExt_ + jy;
        yBufR_[jv] = view(jx, jy, jz);
      }
    }
  }

  // Compute direct transform
  fftw_execute(yPlan_r2c_);

  // Convolution in spectral space
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jl = 0; jl < nl_; ++jl) {
        const size_t jv = jx*(nz_*nl_) + jz*nl_ + jl;
        yBufC_[jv][0] *= ySpecStdDev_[jl];
        yBufC_[jv][1] *= ySpecStdDev_[jl];
      }
    }
  }

  // Compute inverse transform
  fftw_execute(yPlan_c2r_);

  // Deserialize and normalize
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t jy = 0; jy < nyExt_; ++jy) {
        const size_t jv = jx*(nz_*nyExt_) + jz*nyExt_ + jy;
        view(jx, jy, jz) = yBufR_[jv]*yNormFFT_;
      }
    }
  }

  oops::Log::trace() << classname() << "::colsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::vertConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertConvolution starting" << std::endl;

  // Copy field
  atlas::Field copyField = field.clone();
  const auto copyView = atlas::array::make_view<double, 3>(copyField);

  // Apply kernel
  auto view = atlas::array::make_view<double, 3>(field);
  view.assign(0.0);
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jy = 0; jy < nyExt_; ++jy) {
      for (size_t jz = 0; jz < nz_; ++jz) {
        for (size_t jkz = 0; jkz < zKernelSize_; ++jkz) {
          const size_t jzz = jz-jkz+(zKernelSize_-1)/2;
          if (jzz >= 0 && jzz < nz_) {
            view(jx, jy, jz) += copyView(jx, jy, jzz)*zKernel_[jkz];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::vertConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::vertNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertNormalization starting" << std::endl;

  // Apply normalization
  auto view = atlas::array::make_view<double, 3>(field);
  for (size_t jx = 0; jx < nxPerTask_[myrank_]; ++jx) {
    for (size_t jy = 0; jy < nyExt_; ++jy) {
      for (size_t jnz = 0; jnz < zNormSize_; ++jnz) {
        view(jx, jy, jnz) *= zNorm_[jnz];
        view(jx, jy, nz_-1-jnz) *= zNorm_[jnz];
      }
    }
  }

  oops::Log::trace() << classname() << "::vertNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::multiplyRedSqrt(const atlas::Field & colsField,
                                atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrt starting" << std::endl;

  // Create intermediate fields
  atlas::Field colsFieldTmp = colsField.clone();
  atlas::Field rowsField("dummy",
                         atlas::array::make_datatype<double>(),
                         atlas::array::make_shape(nxExt_, nyPerTask_[myrank_], nz_));

  if (nz_ > 1) {
    // Apply vertical kernel
    vertConvolution(colsFieldTmp);

    // Apply vertical normalization
    vertNormalization(colsFieldTmp);
  }

  // Convolution on columns
  colsConvolution(colsFieldTmp);

  // Columns to rows
  colsToRows(colsFieldTmp, rowsField);

  // Apply kernel on rows
  rowsConvolution(rowsField);

  // Rows to reduced grid
  rowsToRed(rowsField, redField);

  oops::Log::trace() << classname() << "::multiplyRedSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerSpec::multiplyRedSqrtTrans(const atlas::Field & redField,
                                     atlas::Field & colsField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans starting" << std::endl;

  // Create intermediate fields
  atlas::Field rowsField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxExt_, nyPerTask_[myrank_], nz_));

  // Reduced grid to rows
  redToRows(redField, rowsField);

  // Convolution on rows
  rowsConvolution(rowsField);

  // Rows to columns
  rowsToCols(rowsField, colsField);

  // Convolution on columns
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
