/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/Layer.h"

#include <netcdf.h>

#include <algorithm>
#include <utility>

#include "atlas/array.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

void Layer::setupVerticalCoord(const atlas::Field & fieldRv,
                               const atlas::Field & fieldWgt) {
  oops::Log::trace() << classname() << "::setupVerticalCoord starting" << std::endl;

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Compute horizontally-averaged vertical length-scale and vertical coordinate
  std::vector<double> vertCoord(nz0_, 0.0);
  std::vector<double> rv(nz0_, 0.0);
  std::vector<double> wgt(nz0_, 0.0);
  const auto viewRv = atlas::array::make_view<double, 2>(fieldRv);
  const auto viewWgt = atlas::array::make_view<double, 2>(fieldWgt);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (viewGhost0(jnode0) == 0) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        double VC = static_cast<double>(k0+1);
        if (gdata_.fieldSet().has_field("vert_coord")) {
          const atlas::Field fieldVertCoord = gdata_.fieldSet()["vert_coord"];
          const auto viewVertCoord = atlas::array::make_view<double, 2>(fieldVertCoord);
          VC = viewVertCoord(jnode0, k0);
        }
        vertCoord[k0] += VC*viewWgt(jnode0, k0);
        rv[k0] += viewRv(jnode0, k0)*viewWgt(jnode0, k0);
        wgt[k0] += viewWgt(jnode0, k0);
      }
    }
  }
  comm_.allReduceInPlace(vertCoord.begin(), vertCoord.end(), eckit::mpi::sum());
  comm_.allReduceInPlace(rv.begin(), rv.end(), eckit::mpi::sum());
  comm_.allReduceInPlace(wgt.begin(), wgt.end(), eckit::mpi::sum());

  // Apply weight
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    ASSERT(wgt[k0] > 0.0);
    vertCoord[k0] = vertCoord[k0]/wgt[k0];
    rv[k0] = rv[k0]/wgt[k0];
  }

  // Compute thickness
  std::vector<double> thickness(nz0_, 0.0);
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    if (k0 == 0) {
      thickness[k0] = std::abs(vertCoord[k0+1]-vertCoord[k0]);
    } else if (k0 == nz0_-1) {
      thickness[k0] = std::abs(vertCoord[k0]-vertCoord[k0-1]);
    } else {
      thickness[k0] = 0.5*std::abs(vertCoord[k0+1]-vertCoord[k0-1]);
    }
  }

  // Normalize thickness with vertical length-scale
  std::vector<double> normThickness(nz0_, 0.0);
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    if (nz0_ > 1) {
      ASSERT(rv[k0] > 0.0);
      normThickness[k0] = thickness[k0]/rv[k0];
    }
  }

  // Compute normalized vertical coordinate
  normVertCoord_.resize(nz0_, 0.0);
  for (size_t k0 = 1; k0 < nz0_; ++k0) {
    normVertCoord_[k0] = normVertCoord_[k0-1]+0.5*(normThickness[k0]+normThickness[k0-1]);
  }

  if (nz0_ > 1) {
    // Rescale normalized vertical coordinate from 0 to nz0_-1
    const double maxNormVertCoord = normVertCoord_[nz0_-1];
    for (size_t k0 = 0; k0 < nz0_; ++k0) {
      normVertCoord_[k0] = normVertCoord_[k0]/maxNormVertCoord*static_cast<double>(nz0_-1);
    }

    // Save rescaled vertical length-scale
    rv_ = static_cast<double>(nz0_-1)/maxNormVertCoord;
  } else {
    // Save rescaled vertical length-scale
    rv_ = 1.0;
  }

  // Copy resolution
  ASSERT(params_.resol.value() != boost::none);
  resol_ = *params_.resol.value();

  oops::Log::trace() << classname() << "::setupVerticalCoord done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::setupInterpolation() {
  oops::Log::trace() << classname() << "::setupInterpolation starting" << std::endl;

  // Model grid indices
  atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
  atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
  auto viewIndexI0 = atlas::array::make_view<int, 1>(fieldIndexI0);
  auto viewIndexJ0 = atlas::array::make_view<int, 1>(fieldIndexJ0);

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Define reduction factors
  double horRedFac = std::max(rh_/resol_, 1.0);
  double verRedFac = std::max(rv_/resol_, 1.0);

  // Reduced grid size
  nx_ = std::min(nx0_, static_cast<size_t>(static_cast<double>(nx0_-1)/horRedFac)+2);
  ny_ = std::min(ny0_, static_cast<size_t>(static_cast<double>(ny0_-1)/horRedFac)+2);
  nz_ = std::min(nz0_, static_cast<size_t>(static_cast<double>(nz0_-1)/verRedFac)+2);
  xRedFac_ = static_cast<double>(nx0_-1)/static_cast<double>(nx_-1);
  yRedFac_ = static_cast<double>(ny0_-1)/static_cast<double>(ny_-1);
  if (nz_ > 1) {
    zRedFac_ = static_cast<double>(nz0_-1)/static_cast<double>(nz_-1);
  } else {
    zRedFac_ = 1.0;
  }

  oops::Log::info() << "Info     :     Target reduction factors: " << std::endl;
  oops::Log::info() << "Info     :     - horizontal: " << horRedFac << std::endl;
  oops::Log::info() << "Info     :     - vertical: " << verRedFac << std::endl;
  oops::Log::info() << "Info     :     Real reduction factors: " << std::endl;
  oops::Log::info() << "Info     :     - along x: " << xRedFac_ << std::endl;
  oops::Log::info() << "Info     :     - along y: " << yRedFac_ << std::endl;
  oops::Log::info() << "Info     :     - along z: " << zRedFac_ << std::endl;

  // Reduced grid coordinates
  std::vector<double> xCoord;
  for (size_t i = 0; i < nx_; ++i) {
    xCoord.push_back(static_cast<double>(i*(nx0_-1))/static_cast<double>(nx_-1));
  }
  std::vector<double> yCoord;
  for (size_t j = 0; j < ny_; ++j) {
    yCoord.push_back(static_cast<double>(j*(ny0_-1))/static_cast<double>(ny_-1));
  }
  std::vector<double> zCoord;
  if (nz_ > 1) {
    for (size_t k = 0; k < nz_; ++k) {
      zCoord.push_back(static_cast<double>(k*(nz0_-1))/static_cast<double>(nz_-1));
    }
  } else {
    zCoord.push_back(0.0);
  }

  // Define reduced grid horizontal distribution
  std::vector<int> mpiTask(nx_*ny_, -1);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (viewGhost0(jnode0) == 0) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t i = 0; i < nx_; ++i) {
          if (viewIndexI0(jnode0)-1 == std::round(xCoord[i]) &&
              viewIndexJ0(jnode0)-1 == std::round(yCoord[j])) {
            mpiTask[i*ny_+j] = myrank_;
          }
        }
      }
    }
  }

  // Check that every point is assigned to a task and remove verification offset
  comm_.allReduceInPlace(mpiTask.begin(), mpiTask.end(), eckit::mpi::sum());
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (mpiTask[i*ny_+j] == -static_cast<int>(comm_.size())) {
        oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")" << std::endl;
        throw eckit::Exception("task not define for this point", Here());
      }
      mpiTask[i*ny_+j] += comm_.size()-1;
    }
  }

  // Create reduced grid FunctionSpace on each task
  std::vector<atlas::PointXY> v;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (static_cast<size_t>(mpiTask[i*ny_+j]) == myrank_) {
        atlas::PointXY p({xCoord[i]/static_cast<double>(nx0_-1),
                          yCoord[j]/static_cast<double>(ny0_-1)});  // Fake coordinate in [0,1]
        v.push_back(p);
      }
    }
  }
  nodes_ = v.size();
  fspace_ = atlas::functionspace::PointCloud(v);

  // Create reduced grid index fields
  atlas::Field fieldIndexI = fspace_.createField<int>(atlas::option::name("index_i"));
  atlas::Field fieldIndexJ = fspace_.createField<int>(atlas::option::name("index_j"));
  auto viewIndexI = atlas::array::make_view<int, 1>(fieldIndexI);
  auto viewIndexJ = atlas::array::make_view<int, 1>(fieldIndexJ);
  size_t jnode = 0;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (static_cast<size_t>(mpiTask[i*ny_+j]) == myrank_) {
        viewIndexI(jnode) = i+1;
        viewIndexJ(jnode) = j+1;
        ++jnode;
      }
    }
  }
  fset_.add(fieldIndexI);
  fset_.add(fieldIndexJ);

  // RecvCounts and received points list
  rRecvCounts_.resize(comm_.size());
  std::fill(rRecvCounts_.begin(), rRecvCounts_.end(), 0);
  std::vector<int> recvPointList;
  size_t zero = 0;
  for (size_t j = 0; j < ny_; ++j) {
    double jMin = yCoord[std::max(j-1, zero)];
    double jMax = yCoord[std::min(j+1, ny_-1)];
    for (size_t i = 0; i < nx_; ++i) {
      double iMin = xCoord[std::max(i-1, zero)];
      double iMax = xCoord[std::min(i+1, nx_-1)];
      bool pointNeeded = false;
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (viewGhost0(jnode0) == 0) {
          if (iMin <= static_cast<double>(viewIndexI0(jnode0)-1) &&
              static_cast<double>(viewIndexI0(jnode0)-1) <= iMax &&
              jMin <= static_cast<double>(viewIndexJ0(jnode0)-1) &&
              static_cast<double>(viewIndexJ0(jnode0)-1) <= jMax) {
            pointNeeded = true;
          }
        }
      }
      if (pointNeeded) {
        ++rRecvCounts_[mpiTask[i*ny_+j]];
        recvPointList.push_back(i*ny_+j);
      }
    }
  }

  // Buffer size
  mSize_ = recvPointList.size();

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
    rSendCounts_.push_back(rRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  rSize_ = 0;
  for (const auto & n : rSendCounts_) rSize_ += n;

  // SendDispls
  rSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    rSendDispls_.push_back(rSendDispls_[jt]+rSendCounts_[jt]);
  }

  // Ordered received points list
  std::vector<size_t> offset(comm_.size(), 0);
  std::vector<int> recvPointListOrdered(mSize_);
  for (size_t jr = 0; jr < mSize_; ++jr) {
    size_t jt = mpiTask[recvPointList[jr]];
    size_t jro = rRecvDispls_[jt]+offset[jt];
    recvPointListOrdered[jro] = recvPointList[jr];
    ++offset[jt];
  }
  std::vector<int> sendPointList(rSize_);
  comm_.allToAllv(recvPointListOrdered.data(), rRecvCounts_.data(), rRecvDispls_.data(),
                  sendPointList.data(), rSendCounts_.data(), rSendDispls_.data());

  // Mapping for sent points
  sendMapping_.resize(rSize_);
  for (size_t js = 0; js < rSize_; ++js) {
    bool found = false;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (static_cast<size_t>(sendPointList[js]) == (viewIndexI(jnode)-1)*ny_+viewIndexJ(jnode)-1) {
        ASSERT(!found);
        sendMapping_[js] = jnode;
        found = true;
      }
    }
    ASSERT(found);
  }

  // Compute horizontal interpolation
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    // Interpolation element default values
    std::string interpType = "n";
    size_t index1 = nx_;
    size_t index2 = ny_;
    std::vector<std::pair<size_t, double>> operations;
    if (viewGhost0(jnode0) == 0) {
      // Model grid indices
      const size_t i0 = viewIndexI0(jnode0)-1;
      const size_t j0 = viewIndexJ0(jnode0)-1;
      const double di = static_cast<double>(i0)/xRedFac_;
      const double dj = static_cast<double>(j0)/yRedFac_;
      const size_t i = static_cast<size_t>(di);
      const size_t j = static_cast<size_t>(dj);
      const bool integerI = (std::abs(static_cast<double>(i)-di) < 1.0e-12);
      const bool integerJ = (std::abs(static_cast<double>(j)-dj) < 1.0e-12);
      const double alphaI = di - static_cast<double>(i);
      const double alphaJ = dj - static_cast<double>(j);
      index1 = i;
      index2 = j;
      if (integerI && integerJ) {
        // Colocated point
        size_t jroIJ = 0;
        bool foundIJ = false;
        size_t jro = 0;
        while ((!foundIJ) && jro < mSize_) {
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          ++jro;
        }
        ASSERT(foundIJ);
        interpType = "c";
        operations.push_back(std::make_pair(jroIJ, 1.0));
      } else if (integerJ) {
        // Linear interpolation along x
        size_t jroIJ = 0;
        size_t jroIpJ = 0;
        bool foundIJ = false;
        bool foundIpJ = false;
        size_t jro = 0;
        while ((!(foundIJ && foundIpJ)) && jro < mSize_) {
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(recvPointListOrdered[jro]) == (i+1)*ny_+j) {
            ASSERT(!foundIpJ);
            jroIpJ = jro;
            foundIpJ = true;
          }
          ++jro;
        }
        ASSERT(foundIJ);
        ASSERT(foundIpJ);
        interpType = "x";
        operations.push_back(std::make_pair(jroIJ, 1.0-alphaI));
        operations.push_back(std::make_pair(jroIpJ, alphaI));
      } else if (integerI) {
        // Linear interpolation along y
        size_t jroIJ = 0;
        size_t jroIJp = 0;
        bool foundIJ = false;
        bool foundIJp = false;
        size_t jro = 0;
        while ((!(foundIJ && foundIJp)) && jro < mSize_) {
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j+1) {
            ASSERT(!foundIJp);
            jroIJp = jro;
            foundIJp = true;
          }
          ++jro;
        }
        ASSERT(foundIJ);
        ASSERT(foundIJp);
        interpType = "y";
        operations.push_back(std::make_pair(jroIJ, 1.0-alphaJ));
        operations.push_back(std::make_pair(jroIJp, alphaJ));
      } else {
        // Bilinear interpolation
        size_t jroIJ = 0;
        size_t jroIpJ = 0;
        size_t jroIJp = 0;
        size_t jroIpJp = 0;
        bool foundIJ = false;
        bool foundIpJ = false;
        bool foundIJp = false;
        bool foundIpJp = false;
        size_t jro = 0;
        while ((!(foundIJ && foundIpJ && foundIJp && foundIpJp)) && jro < mSize_) {
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(recvPointListOrdered[jro]) == (i+1)*ny_+j) {
            ASSERT(!foundIpJ);
            jroIpJ = jro;
            foundIpJ = true;
          }
          if (static_cast<size_t>(recvPointListOrdered[jro]) == i*ny_+j+1) {
            ASSERT(!foundIJp);
            jroIJp = jro;
            foundIJp = true;
          }
          if (static_cast<size_t>(recvPointListOrdered[jro]) == (i+1)*ny_+j+1) {
            ASSERT(!foundIpJp);
            jroIpJp = jro;
            foundIpJp = true;
          }
          ++jro;
        }
        ASSERT(foundIJ);
        ASSERT(foundIpJ);
        ASSERT(foundIJp);
        ASSERT(foundIpJp);
        interpType = "b";
        operations.push_back(std::make_pair(jroIJ, (1.0-alphaI)*(1.0-alphaJ)));
        operations.push_back(std::make_pair(jroIpJ, alphaI*(1.0-alphaJ)));
        operations.push_back(std::make_pair(jroIJp, (1.0-alphaI)*alphaJ));
        operations.push_back(std::make_pair(jroIpJp, alphaI*alphaJ));
      }
    }
    horInterp_.push_back(InterpElement(interpType, index1, index2, operations));
  }

  // Compute vertical interpolation
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    size_t index1 = nz_;
    std::vector<std::pair<size_t, double>> operations;
    if (k0 == 0) {
      // First level
      index1 = 0;
      operations.push_back(std::make_pair(0, 1.0));
    } else if (k0 == nz0_-1) {
      // Last level
      index1 = nz_-1;
      operations.push_back(std::make_pair(nz_-1, 1.0));
    } else {
      // Other levels
      bool found = false;
      size_t k = 0;
      while ((!found) && (k < nz_-1)) {
        if (zCoord[k] <= normVertCoord_[k0] && normVertCoord_[k0] < zCoord[k+1]) {
          index1 = k;
          const double alphaK = normVertCoord_[k0]-zCoord[k];
          if (alphaK > 0.0) {
            // Linear interpolation
            operations.push_back(std::make_pair(k, 1.0-alphaK));
            operations.push_back(std::make_pair(k+1, alphaK));
          } else {
            // Colocated point
            operations.push_back(std::make_pair(k, 1.0));
          }
          found = true;
        }
        ++k;
      }
      ASSERT(found);
    }
    verInterp_.push_back(InterpElement(index1, operations));
  }

  if (!params_.skipTests.value()) {
    // Check interpolation accuracy
    atlas::Field fieldRed = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field fieldModel = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    auto viewRed = atlas::array::make_view<double, 2>(fieldRed);
    auto viewModel = atlas::array::make_view<double, 2>(fieldModel);
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      const double x = static_cast<double>(viewIndexI(jnode)-1)/static_cast<double>(nx_-1);
      const double y = static_cast<double>(viewIndexJ(jnode)-1)/static_cast<double>(ny_-1);
      for (size_t k = 0; k < nz_; ++k) {
        const double z = zCoord[k]/static_cast<double>(nz0_-1);
        viewRed(jnode, k) = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)
          *std::cos(2.0*M_PI*z)+1.0);
      }
    }
    interpolationTL(fieldRed, fieldModel);
    double accuracy = 0.0;
    double maxVal = 0.0;
    double maxRefVal = 0.0;
    std::vector<double> locMax(3, 0.0);
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        const double x = static_cast<double>(viewIndexI0(jnode0)-1)/static_cast<double>(nx0_-1);
        const double y = static_cast<double>(viewIndexJ0(jnode0)-1)/static_cast<double>(ny0_-1);
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          const double z = normVertCoord_[k0]/static_cast<double>(nz0_-1);
          double refVal = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)*std::cos(2.0*M_PI*z)+1.0);
          double diff = std::abs(viewModel(jnode0, k0)-refVal);
          if (diff > accuracy) {
            accuracy = diff;
            maxVal = viewModel(jnode0, k0);
            maxRefVal = refVal;
            locMax = {x, y, z};
          }
        }
      }
    }
    std::vector<double> accuracyPerTask(comm_.size());
    comm_.allGather(accuracy, accuracyPerTask.begin(), accuracyPerTask.end());
    size_t maxTask = std::distance(accuracyPerTask.begin(), std::max_element(
      accuracyPerTask.begin(), accuracyPerTask.end()));
    comm_.broadcast(maxVal, maxTask);
    comm_.broadcast(maxRefVal, maxTask);
    comm_.broadcast(locMax, maxTask);
    oops::Log::info() << std::setprecision(16)<< "Info     :     Interpolation test accuracy: "
                      << accuracy << " (tolerance = " << params_.accuracyTolerance.value()
                      << "), at " << locMax << " : " << maxVal << " != " << maxRefVal << std::endl;
    oops::Log::test() << "    FastLAM interpolation accuracy test";
    if (accuracy < params_.accuracyTolerance.value()) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("interpolation accuracy test failed for block FastLAM", Here());
    }

    // Check interpolation adjoint

    // Create fields
    atlas::Field fieldRedTL = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field fieldModelTL = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field fieldModelAD = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field fieldRedAD = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    auto viewRedTL = atlas::array::make_view<double, 2>(fieldRedTL);
    auto viewModelTL = atlas::array::make_view<double, 2>(fieldModelTL);
    auto viewModelAD = atlas::array::make_view<double, 2>(fieldModelAD);
    auto viewRedAD = atlas::array::make_view<double, 2>(fieldRedAD);

    // Generate random fields
    size_t seed = 7;  // To avoid impact on future random generator calls
    util::NormalDistribution<double> dist(nodes_*nz_+nodes0_*nz0_, 0.0, 1.0, seed);
    size_t jj = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        viewRedTL(jnode, k) = dist[jj];
        ++jj;
      }
    }
    viewModelAD.assign(0.0);
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          viewModelAD(jnode0, k) = dist[jj];
          ++jj;
        }
      }
    }

    // Interpolation TL/AD
    interpolationTL(fieldRedTL, fieldModelTL);
    interpolationAD(fieldModelAD, fieldRedAD);

    // Adjoint test
    double dp1 = 0.0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        dp1 += viewRedTL(jnode, k)*viewRedAD(jnode, k);
      }
    }
    comm_.allReduceInPlace(dp1, eckit::mpi::sum());
    double dp2 = 0.0;
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (viewGhost0(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          dp2 += viewModelTL(jnode0, k0)*viewModelAD(jnode0, k0);
        }
      }
    }
    comm_.allReduceInPlace(dp2, eckit::mpi::sum());

    // Print result
    oops::Log::info() << std::setprecision(16)
                      << "Info     :     FastLAM interpolation adjoint test: y^t (Ax) = "
                      << dp1 << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                      << params_.adjointTolerance.value() << std::endl;
    oops::Log::test() << "    FastLAM interpolation adjoint test";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < params_.adjointTolerance.value()) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("interpolation adjoint test failed for block FastLAM", Here());
    }
  }

  oops::Log::trace() << classname() << "::setupInterpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::setupParallelization() {
  oops::Log::trace() << classname() << "::setupParallelization starting" << std::endl;

  // Get index fields
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto viewIndexI = atlas::array::make_view<int, 1>(fieldIndexI);
  auto viewIndexJ = atlas::array::make_view<int, 1>(fieldIndexJ);

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

  // Reduced grid <=> rows

  // Define destination task
  xTask_.resize(nodes_);
  xOffset_.resize(nodes_);
  std::vector<int> xOffsetPerTask(comm_.size(), 0);
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    bool found = false;
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      if (static_cast<size_t>(viewIndexJ(jnode))-1 >= nyStart_[jt] &&
          static_cast<size_t>(viewIndexJ(jnode))-1 <= nyEnd_[jt]) {
        xTask_[jnode] = jt;
        xOffset_[jnode] = xOffsetPerTask[jt];
        ++xOffsetPerTask[jt];
        ASSERT(!found);
        found = true;
      }
    }
    ASSERT(found);
  }

  // SendCounts
  xSendCounts_.resize(comm_.size());
  std::fill(xSendCounts_.begin(), xSendCounts_.end(), 0);
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    ++xSendCounts_[xTask_[jnode]];
  }

  // SendDispls
  xSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xSendDispls_.push_back(xSendDispls_[jt]+xSendCounts_[jt]);
  }

  // Allgather SendCounts
  eckit::mpi::Buffer<int> xSendCountsBuffer(comm_.size());
  comm_.allGatherv(xSendCounts_.begin(), xSendCounts_.end(), xSendCountsBuffer);
  std::vector<int> xSendCountsGlb_ = std::move(xSendCountsBuffer.buffer);

  // RecvCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xRecvCounts_.push_back(xSendCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  xSize_ = 0;
  for (const auto & n : xRecvCounts_) xSize_ += n;

  // RecvDispls
  xRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    xRecvDispls_.push_back(xRecvDispls_[jt]+xRecvCounts_[jt]);
  }

  // Communicate indices
  std::vector<int> redIndex_i(nodes_);
  std::vector<int> redIndex_j(nodes_);
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    redIndex_i[jnode] = viewIndexI(jnode)-1;
    redIndex_j[jnode] = viewIndexJ(jnode)-1;
  }
  xIndex_i_.resize(xSize_);
  xIndex_j_.resize(xSize_);
  comm_.allToAllv(redIndex_i.data(), xSendCounts_.data(), xSendDispls_.data(),
                  xIndex_i_.data(), xRecvCounts_.data(), xRecvDispls_.data());
  comm_.allToAllv(redIndex_j.data(), xSendCounts_.data(), xSendDispls_.data(),
                  xIndex_j_.data(), xRecvCounts_.data(), xRecvDispls_.data());

  // Rows <=> columns

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

  // SendCounts
  ySendCounts_.resize(comm_.size());
  std::fill(ySendCounts_.begin(), ySendCounts_.end(), 0);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      ++ySendCounts_[yTask_[j][i]];
    }
  }

  // SendDispls
  ySendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    ySendDispls_.push_back(ySendDispls_[jt]+ySendCounts_[jt]);
  }

  // Allgather SendCounts
  eckit::mpi::Buffer<int> ySendCountsBuffer(comm_.size());
  comm_.allGatherv(ySendCounts_.begin(), ySendCounts_.end(), ySendCountsBuffer);
  std::vector<int> ySendCountsGlb_ = std::move(ySendCountsBuffer.buffer);

  // RecvCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    yRecvCounts_.push_back(ySendCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  ySize_ = 0;
  for (const auto & n : yRecvCounts_) ySize_ += n;

  // RecvDispls
  yRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    yRecvDispls_.push_back(yRecvDispls_[jt]+yRecvCounts_[jt]);
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
  ASSERT(xIndex_i.size() == xSize_);
  ASSERT(xIndex_j.size() == xSize_);
  yIndex_i_.resize(ySize_);
  yIndex_j_.resize(ySize_);
  comm_.allToAllv(xIndex_i.data(), ySendCounts_.data(), ySendDispls_.data(),
                  yIndex_i_.data(), yRecvCounts_.data(), yRecvDispls_.data());
  comm_.allToAllv(xIndex_j.data(), ySendCounts_.data(), ySendDispls_.data(),
                  yIndex_j_.data(), yRecvCounts_.data(), yRecvDispls_.data());

  if (!params_.skipTests.value()) {
    // Tests

    // Generate fields
    atlas::Field fieldRed = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field fieldRows = atlas::Field("dummy",
      atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));
    atlas::Field fieldCols = atlas::Field("dummy",
      atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
    auto viewRed = atlas::array::make_view<double, 2>(fieldRed);
    auto viewRows = atlas::array::make_view<double, 3>(fieldRows);
    auto viewCols = atlas::array::make_view<double, 3>(fieldCols);

    // Fill field with global horizontal index
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        viewRed(jnode, k) = static_cast<double>(redIndex_i[jnode]*ny_+redIndex_j[jnode]);
      }
    }

    // Test reduced grid <=> rows
    redToRows(fieldRed, fieldRows);
    rowsToRed(fieldRows, fieldRed);

    // Print result
    oops::Log::test() << "    FastLAM redToRows test";
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (viewRows(i, j, k) != static_cast<double>(i*ny_+(j+nyStart_[myrank_]))) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("redToRows test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        if (viewRed(jnode, k) != static_cast<double>(redIndex_i[jnode]*ny_+redIndex_j[jnode])) {
          oops::Log::test() << " failed" << std::endl;
          throw eckit::Exception("redToRows test failed for block FastLAM", Here());
        }
      }
    }
    oops::Log::test() << " passed" << std::endl;

    // Test Rows <=> columns
    rowsToCols(fieldRows, fieldCols);
    colsToRows(fieldCols, fieldRows);

    // Print result
    oops::Log::test() << "    FastLAM rowsToCols test";
    for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (viewCols(i, j, k) != static_cast<double>((i+nxStart_[myrank_])*ny_+j)) {
            oops::Log::test() << " failed" << std::endl;
            throw eckit::Exception("rowsToCols test failed for block FastLAM", Here());
          }
        }
      }
    }
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          if (viewRows(i, j, k) != static_cast<double>(i*ny_+(j+nyStart_[myrank_]))) {
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

void Layer::setupKernels() {
  oops::Log::trace() << classname() << "::setupKernels starting" << std::endl;

  // Get kernels size
  xKernelSize_ = 2*static_cast<size_t>((0.5*rh_+1.0e-12)/xRedFac_)+1;
  yKernelSize_ = 2*static_cast<size_t>((0.5*rh_+1.0e-12)/yRedFac_)+1;
  zKernelSize_ = 2*static_cast<size_t>((0.5*rv_+1.0e-12)/zRedFac_)+1;
/*
  ASSERT(xKernelSize_ >= 3);
  ASSERT(yKernelSize_ >= 3);
  if (nz_ > 1) {
    ASSERT(zKernelSize_ >= 3);
  }
*/

  // Create kernels
  xKernel_.resize(xKernelSize_);
  yKernel_.resize(yKernelSize_);
  zKernel_.resize(zKernelSize_);
  double xAlpha = xRedFac_/(0.5*rh_);
  double xNorm = 0.0;
  for (size_t jk = 0; jk < xKernelSize_; ++jk) {
    int jkc = jk-(xKernelSize_-1)/2;
    if (jkc < 0) {
      xKernel_[jk] = xAlpha*static_cast<double>(jkc)+1.0;
    } else {
      xKernel_[jk] = -xAlpha*static_cast<double>(jkc)+1.0;
    }
    xNorm += xKernel_[jk]*xKernel_[jk];
  }
  double yAlpha = yRedFac_/(0.5*rh_);
  double yNorm = 0.0;
  for (size_t jk = 0; jk < yKernelSize_; ++jk) {
    int jkc = jk-(yKernelSize_-1)/2;
    if (jkc < 0) {
      yKernel_[jk] = yAlpha*static_cast<double>(jkc)+1.0;
    } else {
      yKernel_[jk] = -yAlpha*static_cast<double>(jkc)+1.0;
    }
    yNorm += yKernel_[jk]*yKernel_[jk];
  }
  double zAlpha = zRedFac_/(0.5*rv_);
  double zNorm = 0.0;
  for (size_t jk = 0; jk < zKernelSize_; ++jk) {
    int jkc = jk-(zKernelSize_-1)/2;
    if (jkc < 0) {
      zKernel_[jk] = zAlpha*static_cast<double>(jkc)+1.0;
    } else {
      zKernel_[jk] = -zAlpha*static_cast<double>(jkc)+1.0;
    }
    zNorm += zKernel_[jk]*zKernel_[jk];
  }

  // Normalize kernels
  xNorm = 1.0/std::sqrt(xNorm);
  for (size_t jk = 0; jk < xKernelSize_; ++jk) {
    xKernel_[jk] *= xNorm;
  }
  yNorm = 1.0/std::sqrt(yNorm);
  for (size_t jk = 0; jk < yKernelSize_; ++jk) {
    yKernel_[jk] *= yNorm;
  }
  zNorm = 1.0/std::sqrt(zNorm);
  for (size_t jk = 0; jk < zKernelSize_; ++jk) {
    zKernel_[jk] *= zNorm;
  }

  oops::Log::trace() << classname() << "::setupKernels done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::setupNormalization() {
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
  atlas::Field fieldNorm = gdata_.functionSpace().createField<double>(
    atlas::option::name(myVar_) | atlas::option::levels(nz0_));
  auto viewNorm = atlas::array::make_view<double, 2>(fieldNorm);
  norm_.add(fieldNorm);

  // Cost-efficient normalization

  // Reduced grid indices
  atlas::Field fieldIndexI = fset_["index_i"];
  atlas::Field fieldIndexJ = fset_["index_j"];
  auto viewIndexI = atlas::array::make_view<int, 1>(fieldIndexI);
  auto viewIndexJ = atlas::array::make_view<int, 1>(fieldIndexJ);

  // One level only
  const size_t nzSave = nz_;
  nz_ = 1;

  // Correlation at one reduced grid cell size distance

  // Horizontal convolution
  atlas::Field fieldRedHor = fspace_.createField<double>(atlas::option::name("dummy") |
    atlas::option::levels(nz_));
  atlas::Field fieldColsHor = atlas::Field("dummy",
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto viewRedHor = atlas::array::make_view<double, 2>(fieldRedHor);
  size_t nxHalf = xNormSize_+1;
  size_t nyHalf = yNormSize_+1;
  std::vector<double> horConv(4*nxHalf*nyHalf, 0.0);
  for (size_t i = 0; i < nxHalf; ++i) {
    for (size_t j = 0; j < nyHalf; ++j) {
      // Setup dirac point
      viewRedHor.assign(0.0);
      for (size_t jnode = 0; jnode < nodes_; ++jnode) {
        if (viewIndexI(jnode)-1 == static_cast<int>(i)
          && viewIndexJ(jnode)-1 == static_cast<int>(j)) {
          viewRedHor(jnode, 0) = 1.0;
        }
      }

      // Apply convolution on reduced grid
      multiplyRedSqrtTrans(fieldRedHor, fieldColsHor);
      multiplyRedSqrt(fieldColsHor, fieldRedHor);

      // Gather horizontal convolution data
      viewRedHor = atlas::array::make_view<double, 2>(fieldRedHor);
      for (size_t jnode = 0; jnode < nodes_; ++jnode) {
        if (viewIndexI(jnode)-1 == static_cast<int>(i)
          && viewIndexJ(jnode)-1 == static_cast<int>(j)) {
          horConv[4*(i*nyHalf+j)+0] = viewRedHor(jnode, 0);
        }
        if (viewIndexI(jnode)-1 == static_cast<int>(i+1)
          && viewIndexJ(jnode)-1 == static_cast<int>(j)) {
          horConv[4*(i*nyHalf+j)+1] = viewRedHor(jnode, 0);
        }
        if (viewIndexI(jnode)-1 == static_cast<int>(i)
          && viewIndexJ(jnode)-1 == static_cast<int>(j+1)) {
          horConv[4*(i*nyHalf+j)+2] = viewRedHor(jnode, 0);
        }
        if (viewIndexI(jnode)-1 == static_cast<int>(i+1)
          && viewIndexJ(jnode)-1 == static_cast<int>(j+1)) {
          horConv[4*(i*nyHalf+j)+3] = viewRedHor(jnode, 0);
        }
      }
    }
  }
  comm_.allReduceInPlace(horConv.begin(), horConv.end(), eckit::mpi::sum());

  // Reset number of levels
  nz_ = nzSave;

  // One grid point only
  const size_t nxSave = nx_;
  const size_t nxPerTaskSave = nxPerTask_[myrank_];
  const size_t nySave = ny_;
  const size_t nyPerTaskSave = nyPerTask_[myrank_];
  nx_ = 1;
  nxPerTask_[myrank_] = 1;
  ny_ = 1;
  nyPerTask_[myrank_] = 1;

  // Vertical convolution
  atlas::Field fieldColsVer = atlas::Field("dummy",
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto viewColsVer = atlas::array::make_view<double, 3>(fieldColsVer);

  std::vector<double> verConv(2*(nz_-1), 0.0);
  for (size_t k = 0; k < nz_-1; ++k) {
    // Setup dirac point
    viewColsVer.assign(0.0);
    viewColsVer(0, 0, k) = 1.0;

    // Apply vertical normalization
    vertNormalization(fieldColsVer);

    // Apply vertical kernel
    vertConvolution(fieldColsVer);
    vertConvolution(fieldColsVer);

    // Apply vertical normalization
    vertNormalization(fieldColsVer);

    // Gather horizontal convolution data
    verConv[2*k+0] = viewColsVer(0, 0, k);
    verConv[2*k+1] = viewColsVer(0, 0, k+1);
  }

  // Reset number of grid points
  nx_ = nxSave;
  nxPerTask_[myrank_] = nxPerTaskSave;
  ny_ = nySave;
  nyPerTask_[myrank_] = nyPerTaskSave;

  // Ghost points
  const auto viewGhost0 = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Compute horizontal normalization
  viewNorm.assign(0.0);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (viewGhost0(jnode0) == 0) {
      // Define offset
      size_t offsetI = std::min(std::min(horInterp_[jnode0].index1(),
        nx_-2-horInterp_[jnode0].index1()), nxHalf-1);
      size_t offsetJ = std::min(std::min(horInterp_[jnode0].index2(),
        ny_-2-horInterp_[jnode0].index2()), nyHalf-1);
      size_t horOffset = 4*(offsetI*nyHalf+offsetJ);
      ASSERT(horOffset <= 4*(nxHalf*nyHalf-1));

      if (horInterp_[jnode0].interpType() == "c") {
        // Colocated point, no normalization needed
        viewNorm(jnode0, 0) = 1.0;
      } else if (horInterp_[jnode0].interpType() == "x") {
        // Linear interpolation along x
        double xW = horConv[horOffset+0]*horInterp_[jnode0].operations()[0].second
                   +horConv[horOffset+1]*horInterp_[jnode0].operations()[1].second;
        double xE = horConv[horOffset+1]*horInterp_[jnode0].operations()[0].second
                   +horConv[horOffset+0]*horInterp_[jnode0].operations()[1].second;
        viewNorm(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xW
                             +horInterp_[jnode0].operations()[1].second*xE;
      } else if (horInterp_[jnode0].interpType() == "y") {
        // Linear interpolation along y
        double xS = horConv[horOffset+0]*horInterp_[jnode0].operations()[0].second
                   +horConv[horOffset+2]*horInterp_[jnode0].operations()[1].second;
        double xN = horConv[horOffset+2]*horInterp_[jnode0].operations()[0].second
                   +horConv[horOffset+0]*horInterp_[jnode0].operations()[1].second;
        viewNorm(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xS
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
        viewNorm(jnode0, 0) = horInterp_[jnode0].operations()[0].second*xSW
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
  std::vector<double> verNorm(nz0_, 0.0);
  verNorm[0] = 1.0;
  verNorm[nz0_-1] = 1.0;
  for (size_t k0 = 1; k0 < nz0_-1; ++k0) {
    size_t verOffset = 2*verInterp_[k0].index1();
    ASSERT(verOffset < 2*(nz_-1));
    double xB = verConv[verOffset+0]*verInterp_[k0].operations()[0].second
               +verConv[verOffset+1]*verInterp_[k0].operations()[1].second;
    double xT = verConv[verOffset+1]*verInterp_[k0].operations()[0].second
               +verConv[verOffset+0]*verInterp_[k0].operations()[1].second;
    verNorm[k0] = verInterp_[k0].operations()[0].second*xB
                 +verInterp_[k0].operations()[1].second*xT;
  }

  // Compute 3D normalization
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (viewGhost0(jnode0) == 0) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        viewNorm(jnode0, k0) = viewNorm(jnode0, 0)*verNorm[k0];
      }
    }
  }

  // Get normalization factor
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    for (size_t k0 = 0; k0 < nz0_; ++k0) {
      if (viewNorm(jnode0, k0) > 0.0) {
        viewNorm(jnode0, k0) = 1.0/std::sqrt(viewNorm(jnode0, k0));
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
    atlas::Field fieldModel = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field fieldRed = fspace_.createField<double>(atlas::option::name("dummy") |
      atlas::option::levels(nz_));
    atlas::Field fieldCols = atlas::Field("dummy",
      atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
    auto viewModel = atlas::array::make_view<double, 2>(fieldModel);
    auto viewCols = atlas::array::make_view<double, 3>(fieldCols);

    // Model grid indices
    atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
    atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
    auto viewIndexI0 = atlas::array::make_view<int, 1>(fieldIndexI0);
    auto viewIndexJ0 = atlas::array::make_view<int, 1>(fieldIndexJ0);

    // Compute normalization accuracy
    oops::Log::info() << "Info     :     Compute exact normalization" << std::endl;
    atlas::Field fieldNormAcc = gdata_.functionSpace().createField<double>(
      atlas::option::name(myVar_) | atlas::option::levels(nz0_));
    auto viewNormAcc = atlas::array::make_view<double, 2>(fieldNormAcc);
    viewNormAcc.assign(util::missingValue<double>());
    double normAccMax = 0.0;;

    for (size_t i0 = 0; i0 < nx0_; i0 += params_.normAccStride.value()) {
      for (size_t j0 = 0; j0 < ny0_; j0 += params_.normAccStride.value()) {
        for (size_t k0 = 0; k0 < nz0_; k0 += params_.normAccStride.value()) {
          // Set Dirac point
          viewModel.assign(0.0);
          int myJnode0 = -1;
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (viewIndexI0(jnode0)-1 == static_cast<int>(i0)
              && viewIndexJ0(jnode0)-1 == static_cast<int>(j0)) {
              viewModel(jnode0, k0) = 1.0;
              myJnode0 = jnode0;
            }
          }

          // Interpolation AD
          interpolationAD(fieldModel, fieldRed);

          // Adjoint square-root multiplication
          multiplyRedSqrtTrans(fieldRed, fieldCols);

          // Compute exact normalization
          double exactNorm = 0.0;
          for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
            for (size_t j = 0; j < ny_; ++j) {
              for (size_t k = 0; k < nz_; ++k) {
                exactNorm += viewCols(i, j, k)*viewCols(i, j, k);
              }
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
              viewNormAcc(myJnode0, k0) = (viewNorm(myJnode0, k0)-exactNorm)/exactNorm;
              normAccMax = std::max(normAccMax, std::abs(viewNormAcc(myJnode0, k0)));
            }
          }
        }
      }
    }
    normAcc_.add(fieldNormAcc);
    comm_.allReduceInPlace(normAccMax, eckit::mpi::max());
    oops::Log::info() << "Info     :     Cost-effective normalization maximum error: "
      << normAccMax << std::endl;
  }

  oops::Log::trace() << classname() << "::setupNormalization done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::read(const int & id) {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // NetCDF ids
  int retval, nz0_id, xKernelSize_id, yKernelSize_id, zKernelSize_id,
    xNormSize_id, yNormSize_id, zNormSize_id,
    normVertCoord_id, xKernel_id, yKernel_id, zKernel_id,
    xNorm_id, yNorm_id, zNorm_id;
  size_t nz0;

  // Get dimensions ids
  if ((retval = nc_inq_dimid(id, "nz0", &nz0_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "xKernelSize", &xKernelSize_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "yKernelSize", &yKernelSize_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "zKernelSize", &zKernelSize_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "xNormSize", &xNormSize_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "yNormSize", &yNormSize_id))) ERR(retval);
  if ((retval = nc_inq_dimid(id, "zNormSize", &zNormSize_id))) ERR(retval);

  // Get dimensions
  if ((retval = nc_inq_dimlen(id, nz0_id, &nz0))) ERR(retval);
  ASSERT(nz0 == nz0_);
  if ((retval = nc_inq_dimlen(id, xKernelSize_id, &xKernelSize_))) ERR(retval);
  if ((retval = nc_inq_dimlen(id, yKernelSize_id, &yKernelSize_))) ERR(retval);
  if ((retval = nc_inq_dimlen(id, zKernelSize_id, &zKernelSize_))) ERR(retval);
  if ((retval = nc_inq_dimlen(id, xNormSize_id, &xNormSize_))) ERR(retval);
  if ((retval = nc_inq_dimlen(id, yNormSize_id, &yNormSize_))) ERR(retval);
  if ((retval = nc_inq_dimlen(id, zNormSize_id, &zNormSize_))) ERR(retval);

  // Get variables ids
  if ((retval = nc_inq_varid(id, "normVertCoord", &normVertCoord_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "xKernel", &xKernel_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "yKernel", &yKernel_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "zKernel", &zKernel_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "xNorm", &xNorm_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "yNorm", &yNorm_id))) ERR(retval);
  if ((retval = nc_inq_varid(id, "zNorm", &zNorm_id))) ERR(retval);

  // Create arrays
  double normVertCoord[nz0_][1];
  double xKernel[xKernelSize_][1];
  double yKernel[yKernelSize_][1];
  double zKernel[zKernelSize_][1];
  double xNorm[xNormSize_][1];
  double yNorm[yNormSize_][1];
  double zNorm[zNormSize_][1];

  // Read data
  if ((retval = nc_get_var_double(id, normVertCoord_id, &normVertCoord[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, xKernel_id, &xKernel[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, yKernel_id, &yKernel[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, zKernel_id, &zKernel[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, xNorm_id, &xNorm[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, yNorm_id, &yNorm[0][0]))) ERR(retval);
  if ((retval = nc_get_var_double(id, zNorm_id, &zNorm[0][0]))) ERR(retval);

  // Resize vectors
  normVertCoord_.resize(nz0_);
  xKernel_.resize(xKernelSize_);
  yKernel_.resize(yKernelSize_);
  zKernel_.resize(zKernelSize_);
  xNorm_.resize(xNormSize_);
  yNorm_.resize(yNormSize_);
  zNorm_.resize(zNormSize_);

  // Copy data
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    normVertCoord_[k0] = normVertCoord[k0][0];
  }
  for (size_t jk = 0; jk < xKernelSize_; ++jk) {
    xKernel_[jk] = xKernel[jk][0];
  }
  for (size_t jk = 0; jk < yKernelSize_; ++jk) {
    yKernel_[jk] = yKernel[jk][0];
  }
  for (size_t jk = 0; jk < zKernelSize_; ++jk) {
    zKernel_[jk] = zKernel[jk][0];
  }
  for (size_t jn = 0; jn < xNormSize_; ++jn) {
    xNorm_[jn] = xNorm[jn][0];
  }
  for (size_t jn = 0; jn < yNormSize_; ++jn) {
    yNorm_[jn] = yNorm[jn][0];
  }
  for (size_t jn = 0; jn < zNormSize_; ++jn) {
    zNorm_[jn] = zNorm[jn][0];
  }

  // Get rh_, rv_ and resol_ as attributes
  if ((retval = nc_get_att_double(id, NC_GLOBAL, "rh", &rh_))) ERR(retval);
  if ((retval = nc_get_att_double(id, NC_GLOBAL, "rv", &rv_))) ERR(retval);
  if ((retval = nc_get_att_double(id, NC_GLOBAL, "resol", &resol_))) ERR(retval);
  if (params_.resol.value() != boost::none) {
    ASSERT(*params_.resol.value() == resol_);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::broadcast() {
  oops::Log::trace() << classname() << "::broadcast starting" << std::endl;

  // Broadcast sizes from root task
  comm_.broadcast(xKernelSize_, 0);
  comm_.broadcast(yKernelSize_, 0);
  comm_.broadcast(zKernelSize_, 0);
  comm_.broadcast(xNormSize_, 0);
  comm_.broadcast(yNormSize_, 0);
  comm_.broadcast(zNormSize_, 0);

  // Resize vectors
  normVertCoord_.resize(nz0_);
  xKernel_.resize(xKernelSize_);
  yKernel_.resize(yKernelSize_);
  zKernel_.resize(zKernelSize_);
  xNorm_.resize(xNormSize_);
  yNorm_.resize(yNormSize_);
  zNorm_.resize(zNormSize_);

  // Broadcast data from root task
  comm_.broadcast(rh_, 0);
  comm_.broadcast(rv_, 0);
  comm_.broadcast(resol_, 0);
  comm_.broadcast(normVertCoord_.begin(), normVertCoord_.end(), 0);
  comm_.broadcast(xKernel_.begin(), xKernel_.end(), 0);
  comm_.broadcast(yKernel_.begin(), yKernel_.end(), 0);
  comm_.broadcast(zKernel_.begin(), zKernel_.end(), 0);
  comm_.broadcast(xNorm_.begin(), xNorm_.end(), 0);
  comm_.broadcast(yNorm_.begin(), yNorm_.end(), 0);
  comm_.broadcast(zNorm_.begin(), zNorm_.end(), 0);

  oops::Log::trace() << classname() << "::broadcast done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<int> Layer::writeDef(const int & id) const {
  oops::Log::trace() << classname() << "::writeDef starting" << std::endl;

  // Vector of ids
  std::vector<int> varIds(8);

  // Copy group id
  varIds[0] = id;

  // NetCDF ids
  int retval, nz0_id[1], xKernelSize_id[1], yKernelSize_id[1], zKernelSize_id[1],
    xNormSize_id[1], yNormSize_id[1], zNormSize_id[1];

  // Create dimensions
  if ((retval = nc_def_dim(id, "nz0", nz0_, &nz0_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "xKernelSize", xKernelSize_, &xKernelSize_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "yKernelSize", yKernelSize_, &yKernelSize_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "zKernelSize", zKernelSize_, &zKernelSize_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "xNormSize", xNormSize_, &xNormSize_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "yNormSize", yNormSize_, &yNormSize_id[0]))) ERR(retval);
  if ((retval = nc_def_dim(id, "zNormSize", zNormSize_, &zNormSize_id[0]))) ERR(retval);

  // Define variables
  if ((retval = nc_def_var(id, "normVertCoord", NC_DOUBLE, 1, nz0_id, &varIds[1]))) ERR(retval);
  if ((retval = nc_def_var(id, "xKernel", NC_DOUBLE, 1, xKernelSize_id, &varIds[2]))) ERR(retval);
  if ((retval = nc_def_var(id, "yKernel", NC_DOUBLE, 1, yKernelSize_id, &varIds[3]))) ERR(retval);
  if ((retval = nc_def_var(id, "zKernel", NC_DOUBLE, 1, zKernelSize_id, &varIds[4]))) ERR(retval);
  if ((retval = nc_def_var(id, "xNorm", NC_DOUBLE, 1, xNormSize_id, &varIds[5]))) ERR(retval);
  if ((retval = nc_def_var(id, "yNorm", NC_DOUBLE, 1, yNormSize_id, &varIds[6]))) ERR(retval);
  if ((retval = nc_def_var(id, "zNorm", NC_DOUBLE, 1, zNormSize_id, &varIds[7]))) ERR(retval);

  // Put rh_, rv_ and resol_ as attributes
  if ((retval = nc_put_att_double(id, NC_GLOBAL, "rh", NC_DOUBLE, 1, &rh_))) ERR(retval);
  if ((retval = nc_put_att_double(id, NC_GLOBAL, "rv", NC_DOUBLE, 1, &rv_))) ERR(retval);
  if ((retval = nc_put_att_double(id, NC_GLOBAL, "resol", NC_DOUBLE, 1, &resol_))) ERR(retval);

  oops::Log::trace() << classname() << "::writeDef done" << std::endl;
  return varIds;
}

// -----------------------------------------------------------------------------

void Layer::writeData(const std::vector<int> & varIds) const {
  oops::Log::trace() << classname() << "::writeData starting" << std::endl;

  // Create arrays
  double normVertCoord[nz0_][1];
  double xKernel[xKernelSize_][1];
  double yKernel[yKernelSize_][1];
  double zKernel[zKernelSize_][1];
  double xNorm[xNormSize_][1];
  double yNorm[yNormSize_][1];
  double zNorm[zNormSize_][1];

  // Copy data
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    normVertCoord[k0][0] = normVertCoord_[k0];
  }
  for (size_t jk = 0; jk < xKernelSize_; ++jk) {
    xKernel[jk][0] = xKernel_[jk];
  }
  for (size_t jk = 0; jk < yKernelSize_; ++jk) {
    yKernel[jk][0] = yKernel_[jk];
  }
  for (size_t jk = 0; jk < zKernelSize_; ++jk) {
    zKernel[jk][0] = zKernel_[jk];
  }
  for (size_t jn = 0; jn < xNormSize_; ++jn) {
    xNorm[jn][0] = xNorm_[jn];
  }
  for (size_t jn = 0; jn < yNormSize_; ++jn) {
    yNorm[jn][0] = yNorm_[jn];
  }
  for (size_t jn = 0; jn < zNormSize_; ++jn) {
    zNorm[jn][0] = zNorm_[jn];
  }

  // Write data
  int retval;
  if ((retval = nc_put_var_double(varIds[0], varIds[1], &normVertCoord[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[2], &xKernel[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[3], &yKernel[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[4], &zKernel[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[5], &xNorm[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[6], &yNorm[0][0]))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[7], &zNorm[0][0]))) ERR(retval);

  oops::Log::trace() << classname() << "::writeData done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::interpolationTL(const atlas::Field & fieldRed, atlas::Field & fieldModel) const {
  oops::Log::trace() << classname() << "::interpolationTL starting" << std::endl;

  // Field on reduced grid
  const auto viewRed = atlas::array::make_view<double, 2>(fieldRed);

  // Scale counts and displs for all levels
  std::vector<int> rSendCounts3D(comm_.size());
  std::vector<int> rSendDispls3D(comm_.size());
  std::vector<int> rRecvCounts3D(comm_.size());
  std::vector<int> rRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rSendCounts3D[jt] = rSendCounts_[jt]*nz_;
    rSendDispls3D[jt] = rSendDispls_[jt]*nz_;
    rRecvCounts3D[jt] = rRecvCounts_[jt]*nz_;
    rRecvDispls3D[jt] = rRecvDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> rVec(rSize_*nz_);
  for (size_t js = 0; js < rSize_; ++js) {
    size_t jnode = sendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < rSize_*nz_);
      rVec[js*nz_+k] = viewRed(jnode, k);
    }
  }

  // Communication
  std::vector<double> mVec(mSize_*nz_);
  comm_.allToAllv(rVec.data(), rSendCounts3D.data(), rSendDispls3D.data(),
                  mVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data());

  // Field on model grid
  auto viewModel = atlas::array::make_view<double, 2>(fieldModel);

  // Interpolation
  viewModel.assign(0.0);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    for (const auto & horOperation : horInterp_[jnode0].operations()) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        for (const auto & verOperation : verInterp_[k0].operations()) {
          size_t mIndex = horOperation.first*nz_+verOperation.first;
          viewModel(jnode0, k0) += horOperation.second*verOperation.second*mVec[mIndex];
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::interpolationTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::interpolationAD(const atlas::Field & fieldModel, atlas::Field & fieldRed) const {
  oops::Log::trace() << classname() << "::interpolationAD starting" << std::endl;

  // Field on model grid
  const auto viewModel = atlas::array::make_view<double, 2>(fieldModel);

  // Scale counts and displs for all levels
  std::vector<int> rSendCounts3D(comm_.size());
  std::vector<int> rSendDispls3D(comm_.size());
  std::vector<int> rRecvCounts3D(comm_.size());
  std::vector<int> rRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rSendCounts3D[jt] = rSendCounts_[jt]*nz_;
    rSendDispls3D[jt] = rSendDispls_[jt]*nz_;
    rRecvCounts3D[jt] = rRecvCounts_[jt]*nz_;
    rRecvDispls3D[jt] = rRecvDispls_[jt]*nz_;
  }

  // Interpolation adjoint
  std::vector<double> mVec(mSize_*nz_, 0.0);
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    for (const auto & horOperation : horInterp_[jnode0].operations()) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        for (const auto & verOperation : verInterp_[k0].operations()) {
          size_t mIndex = horOperation.first*nz_+verOperation.first;
          mVec[mIndex] += horOperation.second*verOperation.second*viewModel(jnode0, k0);
        }
      }
    }
  }

  // Communication
  std::vector<double> rVec(rSize_*nz_);
  comm_.allToAllv(mVec.data(), rRecvCounts3D.data(), rRecvDispls3D.data(),
                  rVec.data(), rSendCounts3D.data(), rSendDispls3D.data());

  // Field on reduced grid
  auto viewRed = atlas::array::make_view<double, 2>(fieldRed);

  // Deserialize
  viewRed.assign(0.0);
  for (size_t js = 0; js < rSize_; ++js) {
    size_t jnode = sendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      viewRed(jnode, k) += rVec[js*nz_+k];
    }
  }

  oops::Log::trace() << classname() << "::interpolationAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::redToRows(const atlas::Field & fieldRed, atlas::Field & fieldRows) const {
  oops::Log::trace() << classname() << "::redToRows starting" << std::endl;

  // Field on reduced grid
  const auto viewRed = atlas::array::make_view<double, 2>(fieldRed);

  // Scale counts and displs for all levels
  std::vector<int> xSendCounts3D(comm_.size());
  std::vector<int> xSendDispls3D(comm_.size());
  std::vector<int> xRecvCounts3D(comm_.size());
  std::vector<int> xRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xSendCounts3D[jt] = xSendCounts_[jt]*nz_;
    xSendDispls3D[jt] = xSendDispls_[jt]*nz_;
    xRecvCounts3D[jt] = xRecvCounts_[jt]*nz_;
    xRecvDispls3D[jt] = xRecvDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> vecRed(nodes_*nz_);
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = xSendDispls3D[xTask_[jnode]] + xOffset_[jnode]*nz_ + k;
      ASSERT(jv >= 0);
      ASSERT(jv < nodes_*nz_);
      vecRed[jv] = viewRed(jnode, k);
    }
  }

  // Communication
  std::vector<double> xField(xSize_*nz_);
  comm_.allToAllv(vecRed.data(), xSendCounts3D.data(), xSendDispls3D.data(),
                  xField.data(), xRecvCounts3D.data(), xRecvDispls3D.data());

  // Field on rows
  auto viewRows = atlas::array::make_view<double, 3>(fieldRows);

  // Deserialize
  for (size_t jx = 0; jx < xSize_; ++jx) {
    ASSERT(xIndex_i_[jx] >= 0);
    ASSERT(xIndex_j_[jx] >= static_cast<int>(nyStart_[myrank_]));
    size_t i = xIndex_i_[jx];
    size_t j = xIndex_j_[jx]-nyStart_[myrank_];
    ASSERT(i < nx_);
    ASSERT(j < nyPerTask_[myrank_]);
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = jx*nz_ + k;
      ASSERT(jv < xSize_*nz_);
      viewRows(i, j, k) = xField[jv];
    }
  }

  oops::Log::trace() << classname() << "::redToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::rowsToRed(const atlas::Field & fieldRows, atlas::Field & fieldRed) const {
  oops::Log::trace() << classname() << "::rowsToRed starting" << std::endl;

  // Field on rows
  const auto viewRows = atlas::array::make_view<double, 3>(fieldRows);

  // Scale counts and displs for all levels
  std::vector<int> xSendCounts3D(comm_.size());
  std::vector<int> xSendDispls3D(comm_.size());
  std::vector<int> xRecvCounts3D(comm_.size());
  std::vector<int> xRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    xSendCounts3D[jt] = xSendCounts_[jt]*nz_;
    xSendDispls3D[jt] = xSendDispls_[jt]*nz_;
    xRecvCounts3D[jt] = xRecvCounts_[jt]*nz_;
    xRecvDispls3D[jt] = xRecvDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> xField(xSize_*nz_);
  for (size_t jx = 0; jx < xSize_; ++jx) {
    ASSERT(xIndex_i_[jx] >= 0);
    ASSERT(xIndex_j_[jx] >= static_cast<int>(nyStart_[myrank_]));
    size_t i = xIndex_i_[jx];
    size_t j = xIndex_j_[jx]-nyStart_[myrank_];
    ASSERT(i < nx_);
    ASSERT(j < nyPerTask_[myrank_]);
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = jx*nz_ + k;
      ASSERT(jv >= 0);
      ASSERT(jv < xSize_*nz_);
      xField[jv] = viewRows(i, j, k);
    }
  }

  // Communication
  std::vector<double> vecRed;
  vecRed.resize(nodes_*nz_);
  comm_.allToAllv(xField.data(), xRecvCounts3D.data(), xRecvDispls3D.data(),
                  vecRed.data(), xSendCounts3D.data(), xSendDispls3D.data());

  // Field on reduced grid
  auto viewRed = atlas::array::make_view<double, 2>(fieldRed);

  // Deserialize
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = xSendDispls3D[xTask_[jnode]] + xOffset_[jnode]*nz_ + k;
      ASSERT(jv < nodes_*nz_);
      viewRed(jnode, k) = vecRed[jv];
    }
  }

  oops::Log::trace() << classname() << "::rowsToRed done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::rowsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsConvolution starting" << std::endl;

  // Field on rows
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nx_);
  ASSERT(static_cast<size_t>(field.shape(1)) == nyPerTask_[myrank_]);

  // Copy field
  atlas::Field fieldCopy = atlas::Field(field.name(),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));
  auto viewCopy = atlas::array::make_view<double, 3>(fieldCopy);
  for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t k = 0; k < nz_; ++k) {
        viewCopy(i, j, k) = view(i, j, k);
      }
    }
  }

  // Apply kernel
  view.assign(0.0);
  for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t jk = 0; jk < xKernelSize_; ++jk) {
        size_t ii = i-jk+(xKernelSize_-1)/2;
        if (ii >= 0 && ii < nx_) {
          for (size_t k = 0; k < nz_; ++k) {
            view(i, j, k) += viewCopy(ii, j, k)*xKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::rowsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::rowsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::rowsNormalization starting" << std::endl;

  // Field on rows
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nx_);
  ASSERT(static_cast<size_t>(field.shape(1)) == nyPerTask_[myrank_]);

  // Apply normalization
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

void Layer::rowsToCols(const atlas::Field & fieldRows, atlas::Field & fieldCols) const {
  oops::Log::trace() << classname() << "::rowsToCols starting" << std::endl;

  // Field on rows
  const auto viewRows = atlas::array::make_view<double, 3>(fieldRows);

  // Scale counts and displs for all levels
  std::vector<int> ySendCounts3D(comm_.size());
  std::vector<int> ySendDispls3D(comm_.size());
  std::vector<int> yRecvCounts3D(comm_.size());
  std::vector<int> yRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ySendCounts3D[jt] = ySendCounts_[jt]*nz_;
    ySendDispls3D[jt] = ySendDispls_[jt]*nz_;
    yRecvCounts3D[jt] = yRecvCounts_[jt]*nz_;
    yRecvDispls3D[jt] = yRecvDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> xField(xSize_*nz_);
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        size_t jv = ySendDispls3D[yTask_[j][i]] + yOffset_[j][i]*nz_ + k;
        ASSERT(jv < xSize_*nz_);
        xField[jv] = viewRows(i, j, k);
      }
    }
  }

  // Communication
  std::vector<double> yField(ySize_*nz_);
  comm_.allToAllv(xField.data(), ySendCounts3D.data(), ySendDispls3D.data(),
                  yField.data(), yRecvCounts3D.data(), yRecvDispls3D.data());

  // Field on columns
  auto viewCols = atlas::array::make_view<double, 3>(fieldCols);

  // Deserialize
  for (size_t jy = 0; jy < ySize_; ++jy) {
    ASSERT(yIndex_i_[jy] >= static_cast<int>(nxStart_[myrank_]));
    ASSERT(yIndex_j_[jy] >= 0);
    size_t i = yIndex_i_[jy]-nxStart_[myrank_];
    size_t j = yIndex_j_[jy];
    ASSERT(i < nxPerTask_[myrank_]);
    ASSERT(j < ny_);
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = jy*nz_ + k;
      ASSERT(jv < ySize_*nz_);
      viewCols(i, j, k) = yField[jv];
    }
  }

  oops::Log::trace() << classname() << "::rowsToCols done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::colsToRows(const atlas::Field & fieldCols, atlas::Field & fieldRows) const {
  oops::Log::trace() << classname() << "::colsToRows starting" << std::endl;

  // Field on columns
  const auto viewCols = atlas::array::make_view<double, 3>(fieldCols);

  // Scale counts and displs for all levels
  std::vector<int> ySendCounts3D(comm_.size());
  std::vector<int> ySendDispls3D(comm_.size());
  std::vector<int> yRecvCounts3D(comm_.size());
  std::vector<int> yRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    ySendCounts3D[jt] = ySendCounts_[jt]*nz_;
    ySendDispls3D[jt] = ySendDispls_[jt]*nz_;
    yRecvCounts3D[jt] = yRecvCounts_[jt]*nz_;
    yRecvDispls3D[jt] = yRecvDispls_[jt]*nz_;
  }

  // Serialize
  std::vector<double> yField(ySize_*nz_);
  for (size_t jy = 0; jy < ySize_; ++jy) {
    ASSERT(yIndex_i_[jy] >= static_cast<int>(nxStart_[myrank_]));
    ASSERT(yIndex_j_[jy] >= 0);
    size_t i = yIndex_i_[jy]-nxStart_[myrank_];
    size_t j = yIndex_j_[jy];
    ASSERT(i < nxPerTask_[myrank_]);
    ASSERT(j < ny_);
    for (size_t k = 0; k < nz_; ++k) {
      size_t jv = jy*nz_ + k;
      ASSERT(jv < ySize_*nz_);
      yField[jv] = viewCols(i, j, k);
    }
  }

  // Communication
  std::vector<double> xField(xSize_*nz_);
  comm_.allToAllv(yField.data(), yRecvCounts3D.data(), yRecvDispls3D.data(),
                  xField.data(), ySendCounts3D.data(), ySendDispls3D.data());

  // Field on rows
  auto viewRows = atlas::array::make_view<double, 3>(fieldRows);

  // Deserialize
  for (size_t i = 0; i < nx_; ++i) {
    for (size_t j = 0; j < nyPerTask_[myrank_]; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        size_t jv = ySendDispls3D[yTask_[j][i]] + yOffset_[j][i]*nz_ + k;
        ASSERT(jv < xSize_*nz_);
        viewRows(i, j, k) = xField[jv];
      }
    }
  }

  oops::Log::trace() << classname() << "::colsToRows done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::colsConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsConvolution starting" << std::endl;

  // Field on columns
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nxPerTask_[myrank_]);
  ASSERT(static_cast<size_t>(field.shape(1)) == ny_);

  // Copy field
  atlas::Field fieldCopy = atlas::Field(field.name(),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto viewCopy = atlas::array::make_view<double, 3>(fieldCopy);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        viewCopy(i, j, k) = view(i, j, k);
      }
    }
  }

  // Apply kernel
  view.assign(0.0);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t jk = 0; jk < yKernelSize_; ++jk) {
        size_t jj = j-jk+(yKernelSize_-1)/2;
        if (jj >= 0 && jj < ny_) {
          for (size_t k = 0; k < nz_; ++k) {
            view(i, j, k) += viewCopy(i, jj, k)*yKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::colsConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::colsNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::colsNormalization starting" << std::endl;

  // Field on columns
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nxPerTask_[myrank_]);
  ASSERT(static_cast<size_t>(field.shape(1)) == ny_);

  // Apply normalization
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

void Layer::vertConvolution(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertConvolution starting" << std::endl;

  // Field on columns
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nxPerTask_[myrank_]);
  ASSERT(static_cast<size_t>(field.shape(1)) == ny_);

  // Copy field
  atlas::Field fieldCopy = atlas::Field(field.name(),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto viewCopy = atlas::array::make_view<double, 3>(fieldCopy);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        viewCopy(i, j, k) = view(i, j, k);
      }
    }
  }

  // Apply kernel
  view.assign(0.0);
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        for (size_t jk = 0; jk < zKernelSize_; ++jk) {
          size_t kk = k-jk+(zKernelSize_-1)/2;
          if (kk >= 0 && kk < nz_) {
            view(i, j, k) += viewCopy(i, j, kk)*zKernel_[jk];
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::vertConvolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::vertNormalization(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::vertNormalization starting" << std::endl;

  // Field on columns
  auto view = atlas::array::make_view<double, 3>(field);

  // Check field shape
  ASSERT(static_cast<size_t>(field.shape(0)) == nxPerTask_[myrank_]);
  ASSERT(static_cast<size_t>(field.shape(1)) == ny_);

  // Apply normalization
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

void Layer::multiplyRedSqrt(const atlas::Field & fieldCols, atlas::Field & fieldRed) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrt starting" << std::endl;

  // Create intermediate fields
  atlas::Field fieldRows = atlas::Field("dummy",
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));
  atlas::Field fieldColsTmp = atlas::Field("dummy",
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));
  auto viewCols = atlas::array::make_view<double, 3>(fieldCols);
  auto viewColsTmp = atlas::array::make_view<double, 3>(fieldColsTmp);

  // Copy input field
  for (size_t i = 0; i < nxPerTask_[myrank_]; ++i) {
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t k = 0; k < nz_; ++k) {
        viewColsTmp(i, j, k) = viewCols(i, j, k);
      }
    }
  }

  if (nz_ > 1) {
    // Apply vertical kernel
    vertConvolution(fieldColsTmp);

    // Apply vertical normalization
    vertNormalization(fieldColsTmp);
  }

  // Apply kernel on columns
  colsConvolution(fieldColsTmp);

  // Apply normalization on columns
  colsNormalization(fieldColsTmp);

  // Columns to rows
  colsToRows(fieldColsTmp, fieldRows);

  // Apply kernel on rows
  rowsConvolution(fieldRows);

  // Apply normalization on rows
  rowsNormalization(fieldRows);

  // Rows to reduced grid
  rowsToRed(fieldRows, fieldRed);

  oops::Log::trace() << classname() << "::multiplyRedSqrt done" << std::endl;
}


// -----------------------------------------------------------------------------

void Layer::multiplyRedSqrtTrans(const atlas::Field & fieldRed, atlas::Field & fieldCols) const {
  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans starting" << std::endl;


  // Create intermediate fields
  atlas::Field fieldRows = atlas::Field("dummy",
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nx_, nyPerTask_[myrank_], nz_));

  // Reduced grid to rows
  redToRows(fieldRed, fieldRows);

  // Apply normalization on rows
  rowsNormalization(fieldRows);

  // Apply kernel on rows
  rowsConvolution(fieldRows);

  // Rows to columns
  rowsToCols(fieldRows, fieldCols);

  // Apply normalization on columns
  colsNormalization(fieldCols);

  // Apply kernel on columns
  colsConvolution(fieldCols);

  if (nz_ > 1) {
    // Apply vertical normalization
    vertNormalization(fieldCols);

    // Apply vertical kernel
    vertConvolution(fieldCols);
  }

  oops::Log::trace() << classname() << "::multiplyRedSqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::multiplySqrt(const atlas::Field & fieldCols, atlas::Field & fieldModel) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Create field on reduced grid
  atlas::Field fieldRed = fspace_.createField<double>(atlas::option::name("dummy") |
                                                      atlas::option::levels(nz_));

  // Square-root multiplication on reduced grid
  multiplyRedSqrt(fieldCols, fieldRed);

  // Interpolation TL
  interpolationTL(fieldRed, fieldModel);

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}


// -----------------------------------------------------------------------------

void Layer::multiplySqrtTrans(const atlas::Field & fieldModel, atlas::Field & fieldCols) const {
  oops::Log::trace() << classname() << "::multiplySqrtTrans starting" << std::endl;

  // Create field on reduced grid
  atlas::Field fieldRed = fspace_.createField<double>(atlas::option::name("dummy") |
                                                      atlas::option::levels(nz_));

  // Interpolation AD
  interpolationAD(fieldModel, fieldRed);

  // Adjoint square-root multiplication on reduced grid
  multiplyRedSqrtTrans(fieldRed, fieldCols);

  oops::Log::trace() << classname() << "::multiplySqrtTrans done" << std::endl;
}

// -----------------------------------------------------------------------------

void Layer::multiply(atlas::Field & fieldModel) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Create field on columns
  atlas::Field fieldCols = atlas::Field(fieldModel.name(),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(nxPerTask_[myrank_], ny_, nz_));

  // Adjoint square-root multiplication
  multiplySqrtTrans(fieldModel, fieldCols);

  // Square-root multiplication
  multiplySqrt(fieldCols, fieldModel);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
