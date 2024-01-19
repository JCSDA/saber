/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/LayerHalo.h"

#include <netcdf.h>

#include <algorithm>
#include <utility>

#include "atlas/array.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

static LayerMaker<LayerHalo> makerHalo_("halo");

// -----------------------------------------------------------------------------

void LayerHalo::setupParallelization() {
  oops::Log::trace() << classname() << "::setupParallelization starting" << std::endl;

  // Get index fields
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);

  // Rows convolution

  // Initialize halo mask
  atlas::Field xcPoints("xcPoints", atlas::array::make_datatype<int>(),
    atlas::array::make_shape(nx_, ny_));
  auto xcPointsView = atlas::array::make_view<int, 2>(xcPoints);
  xcPointsView.assign(-2);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int i = indexIView(jnode)-1;
    int j = indexJView(jnode)-1;
    for (size_t jk = 0; jk < xKernelSize_; ++jk) {
      size_t ii = i-jk+(xKernelSize_-1)/2;
      if (ii >= 0 && ii < nx_) {
        xcPointsView(ii, j) = -1;
      }
    }
  }

  // Set indices
  xcSize_ = 0;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (xcPointsView(i, j) == -1) {
        xcPointsView(i, j) = xcSize_;
        ++xcSize_;
      }
    }
  }

  // RecvCounts and recv points list
  xcRecvCounts_.resize(comm_.size());
  std::fill(xcRecvCounts_.begin(), xcRecvCounts_.end(), 0);
  std::vector<int> xcRecvPointsList;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (xcPointsView(i, j) >= 0) {
        ++xcRecvCounts_[mpiTask_[i*ny_+j]];
        xcRecvPointsList.push_back(i*ny_+j);
      }
    }
  }

  // Buffer size
  xcRecvSize_ = 0;
  for (const auto & n : xcRecvCounts_) xcRecvSize_ += n;
  ASSERT(xcRecvSize_ == xcSize_);

  // RecvDispls
  xcRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xcRecvDispls_.push_back(xcRecvDispls_[jt]+xcRecvCounts_[jt]);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> xcRecvCountsBuffer(comm_.size());
  comm_.allGatherv(xcRecvCounts_.begin(), xcRecvCounts_.end(), xcRecvCountsBuffer);
  std::vector<int> xcRecvCountsGlb_ = std::move(xcRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xcSendCounts_.push_back(xcRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  xcSendSize_ = 0;
  for (const auto & n : xcSendCounts_) xcSendSize_ += n;

  // RecvDispls
  xcSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xcSendDispls_.push_back(xcSendDispls_[jt]+xcSendCounts_[jt]);
  }

  // Ordered received points list
  std::vector<size_t> xcRecvOffset(comm_.size(), 0);
  std::vector<int> xcRecvPointsListOrdered(xcRecvSize_);
  std::vector<int> xcRecvMapping(xcRecvSize_);
  for (size_t jr = 0; jr < xcRecvSize_; ++jr) {
    size_t jt = mpiTask_[xcRecvPointsList[jr]];
    size_t jro = xcRecvDispls_[jt]+xcRecvOffset[jt];
    ASSERT(jro >= 0);
    ASSERT(jro < xcRecvSize_);
    xcRecvPointsListOrdered[jro] = xcRecvPointsList[jr];
    xcRecvMapping[jr] = jro;
    ++xcRecvOffset[jt];
  }
  std::vector<int> xcSendPointsList(xcSendSize_, 0);
  comm_.allToAllv(xcRecvPointsListOrdered.data(), xcRecvCounts_.data(), xcRecvDispls_.data(),
    xcSendPointsList.data(), xcSendCounts_.data(), xcSendDispls_.data());

  // Mapping for sent points
  xcSendMapping_.resize(xcSendSize_);
  for (size_t js = 0; js < xcSendSize_; ++js) {
    bool found = false;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      if (static_cast<size_t>(xcSendPointsList[js]) ==
        (indexIView(jnode)-1)*ny_+indexJView(jnode)-1) {
        ASSERT(!found);
        xcSendMapping_[js] = jnode;
        found = true;
      }
    }
    ASSERT(found);
  }

  // Convolution operations
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int i = indexIView(jnode)-1;
    int j = indexJView(jnode)-1;
    for (size_t jk = 0; jk < xKernelSize_; ++jk) {
      size_t ii = i-jk+(xKernelSize_-1)/2;
      if (ii >= 0 && ii < nx_) {
        xcOperations_.push_back(std::make_tuple(jnode, xcRecvMapping[xcPointsView(ii, j)],
          xKernel_[jk]));
      }
    }
  }

  // Columns convolution

  // Initialize halo mask
  atlas::Field ycPoints("ycPoints", atlas::array::make_datatype<int>(),
    atlas::array::make_shape(nx_, ny_));
  auto ycPointsView = atlas::array::make_view<int, 2>(ycPoints);
  ycPointsView.assign(-2);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int i = indexIView(jnode)-1;
    int j = indexJView(jnode)-1;
    for (size_t jk = 0; jk < yKernelSize_; ++jk) {
      size_t jj = j-jk+(yKernelSize_-1)/2;
      if (jj >= 0 && jj < ny_) {
        ycPointsView(i, jj) = -1;
      }
    }
  }

  // Set indices
  ycSize_ = 0;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (ycPointsView(i, j) == -1) {
        ycPointsView(i, j) = ycSize_;
        ++ycSize_;
      }
    }
  }

  // RecvCounts and recv points list
  ycRecvCounts_.resize(comm_.size());
  std::fill(ycRecvCounts_.begin(), ycRecvCounts_.end(), 0);
  std::vector<int> ycRecvPointsList;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (ycPointsView(i, j) >= 0) {
        ++ycRecvCounts_[mpiTask_[i*ny_+j]];
        ycRecvPointsList.push_back(i*ny_+j);
      }
    }
  }

  // Buffer size
  ycRecvSize_ = 0;
  for (const auto & n : ycRecvCounts_) ycRecvSize_ += n;
  ASSERT(ycRecvSize_ == ycSize_);

  // RecvDispls
  ycRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    ycRecvDispls_.push_back(ycRecvDispls_[jt]+ycRecvCounts_[jt]);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> ycRecvCountsBuffer(comm_.size());
  comm_.allGatherv(ycRecvCounts_.begin(), ycRecvCounts_.end(), ycRecvCountsBuffer);
  std::vector<int> ycRecvCountsGlb_ = std::move(ycRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ycSendCounts_.push_back(ycRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  ycSendSize_ = 0;
  for (const auto & n : ycSendCounts_) ycSendSize_ += n;

  // RecvDispls
  ycSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    ycSendDispls_.push_back(ycSendDispls_[jt]+ycSendCounts_[jt]);
  }

  // Ordered received points list
  std::vector<size_t> ycRecvOffset(comm_.size(), 0);
  std::vector<int> ycRecvPointsListOrdered(ycRecvSize_);
  std::vector<int> ycRecvMapping(ycRecvSize_);
  for (size_t jr = 0; jr < ycRecvSize_; ++jr) {
    size_t jt = mpiTask_[ycRecvPointsList[jr]];
    size_t jro = ycRecvDispls_[jt]+ycRecvOffset[jt];
    ASSERT(jro >= 0);
    ASSERT(jro < ycRecvSize_);
    ycRecvPointsListOrdered[jro] = ycRecvPointsList[jr];
    ycRecvMapping[jr] = jro;
    ++ycRecvOffset[jt];
  }
  std::vector<int> ycSendPointsList(ycSendSize_, 0);
  comm_.allToAllv(ycRecvPointsListOrdered.data(), ycRecvCounts_.data(), ycRecvDispls_.data(),
    ycSendPointsList.data(), ycSendCounts_.data(), ycSendDispls_.data());

  // Mapping for received points
  ycSendMapping_.resize(ycSendSize_);
  for (size_t js = 0; js < ycSendSize_; ++js) {
    bool found = false;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      if (static_cast<size_t>(ycSendPointsList[js]) ==
        (indexIView(jnode)-1)*ny_+indexJView(jnode)-1) {
        ASSERT(!found);
        ycSendMapping_[js] = jnode;
        found = true;
      }
    }
    ASSERT(found);
  }

  // Convolution operations
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int i = indexIView(jnode)-1;
    int j = indexJView(jnode)-1;
    for (size_t jk = 0; jk < yKernelSize_; ++jk) {
      size_t jj = j-jk+(yKernelSize_-1)/2;
      if (jj >= 0 && jj < ny_) {
        ycOperations_.push_back(std::make_tuple(jnode, ycRecvMapping[ycPointsView(i, jj)],
          yKernel_[jk]));
      }
    }
  }

  oops::Log::trace() << classname() << "::setupParallelization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::setupNormalization() {
  oops::Log::trace() << classname() << "::setupNormalization starting" << std::endl;

  // Boundary normalization

  // Create boundary normalization
  xNormSize_ = (xKernelSize_-1)/2;
  yNormSize_ = (yKernelSize_-1)/2;
  zNormSize_ = (zKernelSize_-1)/2;
  xNorm_.resize(xNormSize_);
  yNorm_.resize(yNormSize_);
  zNorm_.resize(zNormSize_);

  // Compute boundary normalization
  std::fill(xNorm_.begin(), xNorm_.end(), 0.0);
  std::fill(yNorm_.begin(), yNorm_.end(), 0.0);
  std::fill(zNorm_.begin(), zNorm_.end(), 0.0);
  for (size_t jn = 0; jn < xNormSize_; ++jn) {
    for (size_t jk = xNormSize_-jn; jk < xKernelSize_; ++jk) {
      xNorm_[jn] += xKernel_[jk]*xKernel_[jk];
    }
    xNorm_[jn] = 1.0/std::sqrt(xNorm_[jn]);
  }
  for (size_t jn = 0; jn < yNormSize_; ++jn) {
    for (size_t jk = yNormSize_-jn; jk < yKernelSize_; ++jk) {
      yNorm_[jn] += yKernel_[jk]*yKernel_[jk];
    }
    yNorm_[jn] = 1.0/std::sqrt(yNorm_[jn]);
  }
  for (size_t jn = 0; jn < zNormSize_; ++jn) {
    for (size_t jk = zNormSize_-jn; jk < zKernelSize_; ++jk) {
      zNorm_[jn] += zKernel_[jk]*zKernel_[jk];
    }
    zNorm_[jn] = 1.0/std::sqrt(zNorm_[jn]);
  }

  // Full normalization
  atlas::Field normField = gdata_.functionSpace().createField<double>(
    atlas::option::name(myVar_) | atlas::option::levels(nz0_));
  auto normView = atlas::array::make_view<double, 2>(normField);
  norm_.add(normField);

  // Cost-efficient normalization

  // Reduced grid indices
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);

  // One level only
  const size_t nzSave = nz_;
  nz_ = 1;

  // Correlation at one reduced grid cell size distance

  // Horizontal convolution
  atlas::Field redFieldHor = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));
  auto redViewHor = atlas::array::make_view<double, 2>(redFieldHor);
  size_t nxHalf = xNormSize_+1;
  size_t nyHalf = yNormSize_+1;
  std::vector<double> horConv(4*nxHalf*nyHalf, 0.0);
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
      multiplyRedSqrtTrans(redFieldHor);
      multiplyRedSqrt(redFieldHor);

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
  size_t rSizeSave = 0;
  rSizeSave = rSize_;
  rSize_ = 1;

  // Vertical convolution
  atlas::Field colsVerField("dummy", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(rSize_, nz_));
  auto colsVerView = atlas::array::make_view<double, 2>(colsVerField);
  std::vector<double> verConv(2*(nz_-1), 0.0);
  for (size_t k = 0; k < nz_-1; ++k) {
    // Setup dirac point

    colsVerView.assign(0.0);
    colsVerView(0, k) = 1.0;

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Apply vertical kernel
    vertConvolution(colsVerField);
    vertConvolution(colsVerField);

    // Apply vertical normalization
    vertNormalization(colsVerField);

    // Gather horizontal convolution data
    verConv[2*k+0] = colsVerView(0, k);
    verConv[2*k+1] = colsVerView(0, k+1);
  }

  // Reset number of grid points
  rSize_ = rSizeSave;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Compute horizontal normalization
  normView.assign(0.0);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      // Define offset
      size_t offsetI = std::min(std::min(horInterp_[jnode0].index1(),
        nx_-2-horInterp_[jnode0].index1()), nxHalf-1);
      size_t offsetJ = std::min(std::min(horInterp_[jnode0].index2(),
        ny_-2-horInterp_[jnode0].index2()), nyHalf-1);
      size_t horOffset = 4*(offsetI*nyHalf+offsetJ);
      ASSERT(horOffset <= 4*(nxHalf*nyHalf-1));

      if (horInterp_[jnode0].interpType() == "c") {
        // Colocated point, no normalization needed
        normView(jnode0, 0) = 1.0;
      } else if (horInterp_[jnode0].interpType() == "x") {
        // Linear interpolation along x
        double xW = horConv[horOffset+0]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+1]*horInterp_[jnode0].operations()[1].second;
        double xE = horConv[horOffset+1]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+0]*horInterp_[jnode0].operations()[1].second;
        normView(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xW
          +horInterp_[jnode0].operations()[1].second*xE;
      } else if (horInterp_[jnode0].interpType() == "y") {
        // Linear interpolation along y
        double xS = horConv[horOffset+0]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+2]*horInterp_[jnode0].operations()[1].second;
        double xN = horConv[horOffset+2]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+0]*horInterp_[jnode0].operations()[1].second;
        normView(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xS
          +horInterp_[jnode0].operations()[1].second*xN;
      } else if (horInterp_[jnode0].interpType() == "b") {
        // Bilinear interpolation
        double xSW = horConv[horOffset+0]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+1]*horInterp_[jnode0].operations()[1].second
          +horConv[horOffset+2]*horInterp_[jnode0].operations()[2].second
          +horConv[horOffset+3]*horInterp_[jnode0].operations()[3].second;
        double xSE = horConv[horOffset+1]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+0]*horInterp_[jnode0].operations()[1].second
          +horConv[horOffset+3]*horInterp_[jnode0].operations()[2].second
          +horConv[horOffset+2]*horInterp_[jnode0].operations()[3].second;
        double xNW = horConv[horOffset+2]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+3]*horInterp_[jnode0].operations()[1].second
          +horConv[horOffset+0]*horInterp_[jnode0].operations()[2].second
          +horConv[horOffset+1]*horInterp_[jnode0].operations()[3].second;
        double xNE = horConv[horOffset+3]*horInterp_[jnode0].operations()[0].second
          +horConv[horOffset+2]*horInterp_[jnode0].operations()[1].second
          +horConv[horOffset+1]*horInterp_[jnode0].operations()[2].second
          +horConv[horOffset+0]*horInterp_[jnode0].operations()[3].second;
        normView(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xSW
          +horInterp_[jnode0].operations()[1].second*xSE
          +horInterp_[jnode0].operations()[2].second*xNW
          +horInterp_[jnode0].operations()[3].second*xNE;
      } else {
        throw eckit::Exception("wrong interpolation type: " + horInterp_[jnode0].interpType(),
          Here());
      }
    }
  }

  // Compute vertical normalization
  std::vector<double> verNorm(nz0_, 1.0);
  if (nz0_ > 1) {
    verNorm[0] = 1.0;
    verNorm[nz0_-1] = 1.0;
    for (size_t k0 = 1; k0 < nz0_-1; ++k0) {
      if (verInterp_[k0].operations().size() > 1) {
        size_t verOffset = 2*verInterp_[k0].index1();
        ASSERT(verOffset < 2*(nz_-1));
        double xB = verConv[verOffset+0]*verInterp_[k0].operations()[0].second
          +verConv[verOffset+1]*verInterp_[k0].operations()[1].second;
        double xT = verConv[verOffset+1]*verInterp_[k0].operations()[0].second
          +verConv[verOffset+0]*verInterp_[k0].operations()[1].second;
        verNorm[k0] = verInterp_[k0].operations()[0].second*xB
          +verInterp_[k0].operations()[1].second*xT;
      }
    }
  }

  // Compute 3D normalization
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        normView(jnode0, k0) = normView(jnode0, 0)*verNorm[k0];
      }
    }
  }

  // Get normalization factor
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    for (size_t k0 = 0; k0 < nz0_; ++k0) {
      if (normView(jnode0, k0) > 0.0) {
        normView(jnode0, k0) = 1.0/std::sqrt(normView(jnode0, k0));
      }
    }
  }

  // Check whether normalization accuracy should be computed
  bool computeNormAcc = false;
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = params_.outputModelFilesConf.value().get_value_or({});
  for (const auto & conf : outputModelFilesConf) {
    const std::string param = conf.getString("parameter");
    if (param == "normalization accuracy") {
      computeNormAcc = true;
    }
  }

  if (computeNormAcc) {
    // Brute-force normalization to assess cost-effective normalization accuracy

    // Create fields
    atlas::Field modelField = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy") |
      atlas::option::levels(nz_));
    auto modelView = atlas::array::make_view<double, 2>(modelField);

    // Model grid indices
    atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
    atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
    auto indexIView0 = atlas::array::make_view<int, 1>(fieldIndexI0);
    auto indexJView0 = atlas::array::make_view<int, 1>(fieldIndexJ0);

    // Compute normalization accuracy
    oops::Log::info() << "Info     :     Compute exact normalization" << std::endl;
    atlas::Field normAccField = gdata_.functionSpace().createField<double>(
      atlas::option::name(myVar_) | atlas::option::levels(nz0_));
    auto normAccView = atlas::array::make_view<double, 2>(normAccField);
    normAccView.assign(util::missingValue<double>());
    double normAccMax = 0.0;;

    for (size_t i0 = 0; i0 < nx0_; i0 += params_.normAccStride.value()) {
      for (size_t j0 = 0; j0 < ny0_; j0 += params_.normAccStride.value()) {
        for (size_t k0 = 0; k0 < nz0_; k0 += params_.normAccStride.value()) {
          // Set Dirac point
          modelView.assign(0.0);
          int myJnode0 = -1;
          for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
            if (indexIView0(jnode0)-1 == static_cast<int>(i0)
              && indexJView0(jnode0)-1 == static_cast<int>(j0)) {
              modelView(jnode0, k0) = 1.0;
              myJnode0 = jnode0;
            }
          }

          // Interpolation AD
          interpolationAD(modelField, redField);

          // Adjoint square-root multiplication
          multiplyRedSqrtTrans(redField);

          // Compute exact normalization
          double exactNorm = 0.0;
          const auto redView = atlas::array::make_view<double, 2>(redField);
          for (size_t jnode = 0; jnode < rSize_; ++jnode) {
            for (size_t k = 0; k < nz_; ++k) {
              exactNorm += redView(jnode, k)*redView(jnode, k);
            }
          }
          comm_.allReduceInPlace(exactNorm, eckit::mpi::sum());

          // Get exact normalization factor
          if (exactNorm > 0.0) {
            exactNorm = 1.0/std::sqrt(exactNorm);
          }

          // Assess cost-efficient normalization quality
          if (myJnode0 > -1) {
            if (exactNorm > 0.0) {
              normAccView(myJnode0, k0) = (normView(myJnode0, k0)-exactNorm)/exactNorm;
              normAccMax = std::max(normAccMax, std::abs(normAccView(myJnode0, k0)));
            }
          }
        }
      }
    }
    normAcc_.add(normAccField);
    comm_.allReduceInPlace(normAccMax, eckit::mpi::max());
    oops::Log::info() << "Info     :     Cost-effective normalization maximum error: "
      << normAccMax << std::endl;
  }

  oops::Log::trace() << classname() << "::setupNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::multiplySqrt(const atlas::Field & cv,
                             atlas::Field & modelField,
                             const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Create field on reduced grid
  atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy") |
                                                      atlas::option::levels(nz_));

  // Control vector to reduced grid
  const auto cvView = atlas::array::make_view<double, 1>(cv);
  auto redView = atlas::array::make_view<double, 2>(redField);
  size_t index = offset;
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      redView(jnode, k) = cvView(index);
      ++index;
    }
  }

  // Square-root multiplication on reduced grid
  multiplyRedSqrt(redField);

  // Interpolation TL
  interpolationTL(redField, modelField);

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::multiplySqrtTrans(const atlas::Field & modelField,
                                  atlas::Field & cv,
                                  const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtTrans starting" << std::endl;

  // Create field on reduced grid
  atlas::Field redField = fspace_.createField<double>(atlas::option::name("dummy") |
                                                      atlas::option::levels(nz_));

  // Interpolation AD
  interpolationAD(modelField, redField);

  // Adjoint square-root multiplication on reduced grid
  multiplyRedSqrtTrans(redField);

  // Reduced grid to control vector
  const auto redView = atlas::array::make_view<double, 2>(redField);
  auto cvView = atlas::array::make_view<double, 1>(cv);
  size_t index = offset;
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      cvView(index) = redView(jnode, k);
      ++index;
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

void LayerHalo::rowsConvolutionTL(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsConvolutionTL starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> xcSendCounts3D(comm_.size());
  std::vector<int> xcSendDispls3D(comm_.size());
  std::vector<int> xcRecvCounts3D(comm_.size());
  std::vector<int> xcRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xcSendCounts3D[jt] = xcSendCounts_[jt]*nz_;
    xcSendDispls3D[jt] = xcSendDispls_[jt]*nz_;
    xcRecvCounts3D[jt] = xcRecvCounts_[jt]*nz_;
    xcRecvDispls3D[jt] = xcRecvDispls_[jt]*nz_;
  }

  // Serialize
  auto view = atlas::array::make_view<double, 2>(field);
  std::vector<double> xcSendVec(xcSendSize_*nz_);
  for (size_t js = 0; js < xcSendSize_; ++js) {
    size_t jnode = xcSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < xcSendSize_*nz_);
      xcSendVec[js*nz_+k] = view(jnode, k);
    }
  }

  // Communication
  std::vector<double> xcRecvVec(xcSize_*nz_);
  comm_.allToAllv(xcSendVec.data(), xcSendCounts3D.data(), xcSendDispls3D.data(),
    xcRecvVec.data(), xcRecvCounts3D.data(), xcRecvDispls3D.data());

  // Convolution
  view.assign(0.0);
  for (const auto & op : xcOperations_) {
    for (size_t k = 0; k < nz_; ++k) {
      view(std::get<0>(op), k) += xcRecvVec[std::get<1>(op)*nz_+k]*std::get<2>(op);
    }
  }

  oops::Log::trace() << classname() << "::rowsConvolutionTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::rowsConvolutionAD(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsConvolutionAD starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> xcSendCounts3D(comm_.size());
  std::vector<int> xcSendDispls3D(comm_.size());
  std::vector<int> xcRecvCounts3D(comm_.size());
  std::vector<int> xcRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xcSendCounts3D[jt] = xcSendCounts_[jt]*nz_;
    xcSendDispls3D[jt] = xcSendDispls_[jt]*nz_;
    xcRecvCounts3D[jt] = xcRecvCounts_[jt]*nz_;
    xcRecvDispls3D[jt] = xcRecvDispls_[jt]*nz_;
  }

  // Convolution
  auto view = atlas::array::make_view<double, 2>(field);
  std::vector<double> xcRecvVec(xcSize_*nz_, 0.0);
  for (const auto & op : xcOperations_) {
    for (size_t k = 0; k < nz_; ++k) {
      xcRecvVec[std::get<1>(op)*nz_+k] += view(std::get<0>(op), k)*std::get<2>(op);
    }
  }

  // Communication
  std::vector<double> xcSendVec(xcSendSize_*nz_);
  comm_.allToAllv(xcRecvVec.data(), xcRecvCounts3D.data(), xcRecvDispls3D.data(),
    xcSendVec.data(), xcSendCounts3D.data(), xcSendDispls3D.data());

  // Deserialize
  view.assign(0.0);
  for (size_t js = 0; js < xcSendSize_; ++js) {
    size_t jnode = xcSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < xcSendSize_*nz_);
      view(jnode, k) += xcSendVec[js*nz_+k];
    }
  }

  oops::Log::trace() << classname() << "::rowsConvolutionAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::rowsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsNormalization starting" << std::endl;

  // Get index fields
  atlas::Field fieldIndexI = fset_["index_i"];
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);

  // Apply normalization
  auto view = atlas::array::make_view<double, 2>(field);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int i = indexIView(jnode)-1;
    if (i < static_cast<int>(xNormSize_)) {
      for (size_t k = 0; k < nz_; ++k) {
        view(jnode, k) *= xNorm_[i];
      }
    } else if (i >= static_cast<int>(nx_-xNormSize_)) {
      for (size_t k = 0; k < nz_; ++k) {
        view(jnode, k) *= xNorm_[nx_-1-i];
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::colsConvolutionTL(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsConvolutionTL starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> ycSendCounts3D(comm_.size());
  std::vector<int> ycSendDispls3D(comm_.size());
  std::vector<int> ycRecvCounts3D(comm_.size());
  std::vector<int> ycRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ycSendCounts3D[jt] = ycSendCounts_[jt]*nz_;
    ycSendDispls3D[jt] = ycSendDispls_[jt]*nz_;
    ycRecvCounts3D[jt] = ycRecvCounts_[jt]*nz_;
    ycRecvDispls3D[jt] = ycRecvDispls_[jt]*nz_;
  }

  // Serialize
  auto view = atlas::array::make_view<double, 2>(field);
  std::vector<double> ycSendVec(ycSendSize_*nz_);
  for (size_t js = 0; js < ycSendSize_; ++js) {
    size_t jnode = ycSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < ycSendSize_*nz_);
      ycSendVec[js*nz_+k] = view(jnode, k);
    }
  }

  // Communication
  std::vector<double> ycRecvVec(ycSize_*nz_);
  comm_.allToAllv(ycSendVec.data(), ycSendCounts3D.data(), ycSendDispls3D.data(),
    ycRecvVec.data(), ycRecvCounts3D.data(), ycRecvDispls3D.data());

  // Convolution
  view.assign(0.0);
  for (const auto & op : ycOperations_) {
    for (size_t k = 0; k < nz_; ++k) {
      view(std::get<0>(op), k) += ycRecvVec[std::get<1>(op)*nz_+k]*std::get<2>(op);
    }
  }

  oops::Log::trace() << classname() << "::colsConvolutionTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::colsConvolutionAD(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsConvolutionAD starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> ycSendCounts3D(comm_.size());
  std::vector<int> ycSendDispls3D(comm_.size());
  std::vector<int> ycRecvCounts3D(comm_.size());
  std::vector<int> ycRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ycSendCounts3D[jt] = ycSendCounts_[jt]*nz_;
    ycSendDispls3D[jt] = ycSendDispls_[jt]*nz_;
    ycRecvCounts3D[jt] = ycRecvCounts_[jt]*nz_;
    ycRecvDispls3D[jt] = ycRecvDispls_[jt]*nz_;
  }

  // Convolution
  auto view = atlas::array::make_view<double, 2>(field);
  std::vector<double> ycRecvVec(ycSize_*nz_, 0.0);
  for (const auto & op : ycOperations_) {
    for (size_t k = 0; k < nz_; ++k) {
      ycRecvVec[std::get<1>(op)*nz_+k] += view(std::get<0>(op), k)*std::get<2>(op);
    }
  }

  // Communication
  std::vector<double> ycSendVec(ycSendSize_*nz_);
  comm_.allToAllv(ycRecvVec.data(), ycRecvCounts3D.data(), ycRecvDispls3D.data(),
    ycSendVec.data(), ycSendCounts3D.data(), ycSendDispls3D.data());


  // Deserialize
  view.assign(0.0);
  for (size_t js = 0; js < ycSendSize_; ++js) {
    size_t jnode = ycSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < ycSendSize_*nz_);
      view(jnode, k) += ycSendVec[js*nz_+k];
    }
  }

  oops::Log::trace() << classname() << "::colsConvolutionAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::colsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsNormalization starting" << std::endl;

  // Get index fields
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);

  // Apply normalization
  auto view = atlas::array::make_view<double, 2>(field);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    int j = indexJView(jnode)-1;
    if (j < static_cast<int>(yNormSize_)) {
      for (size_t k = 0; k < nz_; ++k) {
        view(jnode, k) *= yNorm_[j];
      }
    } else if (j >= static_cast<int>(ny_-yNormSize_)) {
      for (size_t k = 0; k < nz_; ++k) {
        view(jnode, k) *= yNorm_[ny_-1-j];
      }
    }
  }

  oops::Log::trace() << classname() << "::colsNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::vertConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertConvolution starting" << std::endl;

  // Copy field
  atlas::Field copyField = field.clone();
  const auto copyView = atlas::array::make_view<double, 2>(copyField);

  // Apply kernel
  auto view = atlas::array::make_view<double, 2>(field);
  view.assign(0.0);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      for (size_t jk = 0; jk < zKernelSize_; ++jk) {
        size_t kk = k-jk+(zKernelSize_-1)/2;
        if (kk >= 0 && kk < nz_) {
          view(jnode, k) += copyView(jnode, kk)*zKernel_[jk];
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::vertConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::vertNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertNormalization starting" << std::endl;

  // Apply normalization
  auto view = atlas::array::make_view<double, 2>(field);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t k = 0; k < zNormSize_; ++k) {
      view(jnode, k) *= zNorm_[k];
      view(jnode, nz_-1-k) *= zNorm_[k];
    }
  }

  oops::Log::trace() << classname() << "::vertNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::multiplyRedSqrt(atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrt starting" << std::endl;

  if (nz_ > 1) {
    // Apply vertical kernel
    vertConvolution(redField);

    // Apply vertical normalization
    vertNormalization(redField);
  }

  // Apply kernel on columns
  colsConvolutionTL(redField);

  // Apply normalization on columns
  colsNormalization(redField);

  // Apply kernel on rows
  rowsConvolutionTL(redField);

  // Apply normalization on rows
  rowsNormalization(redField);

  oops::Log::trace() << classname() << "::multiplyRedSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerHalo::multiplyRedSqrtTrans(atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans starting" << std::endl;

  // Apply normalization on rows
  rowsNormalization(redField);

  // Apply kernel on rows
  rowsConvolutionAD(redField);

  // Apply normalization on columns
  colsNormalization(redField);

  // Apply kernel on columns
  colsConvolutionAD(redField);

  if (nz_ > 1) {
    // Apply vertical normalization
    vertNormalization(redField);

    // Apply vertical kernel
    vertConvolution(redField);
  }

  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
