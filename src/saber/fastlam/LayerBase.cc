/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/LayerBase.h"

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

LayerFactory::LayerFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::LayerFactory."
      << std::endl;
    throw eckit::Exception("Element already registered in saber::LayerFactory.",
      Here());
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<LayerBase> LayerFactory::create(
  const FastLAMParametersBase & params,
  const oops::GeometryData & gdata,
  const std::string & myVar,
  const size_t & nx0,
  const size_t & ny0,
  const size_t & nz0) {
  oops::Log::trace() << "LayerBase::create starting" << std::endl;
  const std::string id = params.parallelization.value();
  typename std::map<std::string, LayerFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::LayerFactory." << std::endl;
    throw eckit::UserError("Element does not exist in saber::LayerFactory.", Here());
  }
  std::unique_ptr<LayerBase> ptr =
    jsb->second->make(params, gdata, myVar, nx0, ny0, nz0);
  oops::Log::trace() << "LayerBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

void LayerBase::setupVerticalCoord(const atlas::Field & rvField,
                                   const atlas::Field & wgtField) {
  oops::Log::trace() << classname() << "::setupVerticalCoord starting" << std::endl;

  if (nz0_ == 1) {
    // Compute normalized vertical coordinate
    normVertCoord_.resize(nz0_, 0.0);

    // Save rescaled vertical length-scale
    rv_ = 1.0;
  } else {
    // Ghost points
    const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

    // Compute horizontally-averaged vertical length-scale and vertical coordinate
    std::vector<double> vertCoord(nz0_, 0.0);
    std::vector<double> rv(nz0_, 0.0);
    std::vector<double> wgt(nz0_, 0.0);
    const auto rvView = atlas::array::make_view<double, 2>(rvField);
    const auto wgtView = atlas::array::make_view<double, 2>(wgtField);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          double VC = static_cast<double>(k0+1);
          if (gdata_.fieldSet().has("vert_coord")) {
            const atlas::Field vertCoordField = gdata_.fieldSet()["vert_coord"];
            const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoordField);
            VC = vertCoordView(jnode0, k0);
          }
          vertCoord[k0] += VC*wgtView(jnode0, k0);
          rv[k0] += rvView(jnode0, k0)*wgtView(jnode0, k0);
          wgt[k0] += wgtView(jnode0, k0);
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
      ASSERT(rv[k0] > 0.0);
      normThickness[k0] = thickness[k0]/rv[k0];
    }

    // Compute normalized vertical coordinate
    normVertCoord_.resize(nz0_, 0.0);
    for (size_t k0 = 1; k0 < nz0_; ++k0) {
      normVertCoord_[k0] = normVertCoord_[k0-1]+0.5*(normThickness[k0]+normThickness[k0-1]);
    }

    // Rescale normalized vertical coordinate from 0 to nz0_-1
    const double maxNormVertCoord = normVertCoord_[nz0_-1];
    for (size_t k0 = 0; k0 < nz0_; ++k0) {
      normVertCoord_[k0] = normVertCoord_[k0]/maxNormVertCoord*static_cast<double>(nz0_-1);
    }

    // Save rescaled vertical length-scale
    rv_ = static_cast<double>(nz0_-1)/maxNormVertCoord;
  }

  oops::Log::trace() << classname() << "::setupVerticalCoord done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerBase::setupInterpolation() {
  oops::Log::trace() << classname() << "::setupInterpolation starting" << std::endl;

  // Model grid indices
  atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
  atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
  auto indexIView0 = atlas::array::make_view<int, 1>(fieldIndexI0);
  auto indexJView0 = atlas::array::make_view<int, 1>(fieldIndexJ0);

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Reduced grid size
  nx_ = std::min(nx0_, static_cast<size_t>(static_cast<double>(nx0_-1)/rfh_)+2);
  ny_ = std::min(ny0_, static_cast<size_t>(static_cast<double>(ny0_-1)/rfh_)+2);
  nz_ = std::min(nz0_, static_cast<size_t>(static_cast<double>(nz0_-1)/rfv_)+2);
  xRedFac_ = static_cast<double>(nx0_-1)/static_cast<double>(nx_-1);
  yRedFac_ = static_cast<double>(ny0_-1)/static_cast<double>(ny_-1);
  if (nz_ > 1) {
    zRedFac_ = static_cast<double>(nz0_-1)/static_cast<double>(nz_-1);
  } else {
    zRedFac_ = 1.0;
  }

  oops::Log::info() << "Info     :     Target reduction factors: " << std::endl;
  oops::Log::info() << "Info     :     - horizontal: " << rfh_ << std::endl;
  oops::Log::info() << "Info     :     - vertical: " << rfv_ << std::endl;
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
  mpiTask_.resize(nx_*ny_, 0);
  std::vector<int> mpiMask(nx_*ny_, 0);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t i = 0; i < nx_; ++i) {
          if (indexIView0(jnode0)-1 == std::round(xCoord[i]) &&
              indexJView0(jnode0)-1 == std::round(yCoord[j])) {
            mpiTask_[i*ny_+j] = myrank_;
            mpiMask[i*ny_+j] = 1;
          }
        }
      }
    }
  }

  // Check that every point is assigned to a task
  comm_.allReduceInPlace(mpiTask_.begin(), mpiTask_.end(), eckit::mpi::sum());
  comm_.allReduceInPlace(mpiMask.begin(), mpiMask.end(), eckit::mpi::sum());
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (mpiMask[i*ny_+j] == 0) {
        oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")" << std::endl;
        throw eckit::Exception("task not define for this point", Here());
      }
      if (mpiMask[i*ny_+j] > 1) {
        oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")" << std::endl;
        throw eckit::Exception("task defined more than once for this point", Here());
      }
    }
  }

  // Create reduced grid FunctionSpace on each task
  std::vector<atlas::PointXY> v;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (static_cast<size_t>(mpiTask_[i*ny_+j]) == myrank_) {
        atlas::PointXY p({xCoord[i]/static_cast<double>(nx0_-1),
                          yCoord[j]/static_cast<double>(ny0_-1)});  // Fake coordinate in [0,1]
        v.push_back(p);
      }
    }
  }
  rSize_ = v.size();
  fspace_ = atlas::functionspace::PointCloud(v);

  // Create reduced grid index fields
  atlas::Field fieldIndexI = fspace_.createField<int>(atlas::option::name("index_i"));
  atlas::Field fieldIndexJ = fspace_.createField<int>(atlas::option::name("index_j"));
  auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
  auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);
  size_t jnode = 0;
  for (size_t j = 0; j < ny_; ++j) {
    for (size_t i = 0; i < nx_; ++i) {
      if (static_cast<size_t>(mpiTask_[i*ny_+j]) == myrank_) {
        indexIView(jnode) = i+1;
        indexJView(jnode) = j+1;
        ++jnode;
      }
    }
  }
  fset_.add(fieldIndexI);
  fset_.add(fieldIndexJ);

  // RecvCounts and received points list
  mRecvCounts_.resize(comm_.size());
  std::fill(mRecvCounts_.begin(), mRecvCounts_.end(), 0);
  std::vector<int> mRecvPointsList;
  for (size_t j = 0; j < ny_; ++j) {
    double jMin = j > 0 ? yCoord[j-1] : yCoord[0];
    double jMax = yCoord[std::min(j+1, ny_-1)];
    for (size_t i = 0; i < nx_; ++i) {
      double iMin = i > 0 ? xCoord[i-1] : xCoord[0];
      double iMax = xCoord[std::min(i+1, nx_-1)];
      bool pointsNeeded = false;
      for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
        if (ghostView(jnode0) == 0) {
          if (iMin <= static_cast<double>(indexIView0(jnode0)-1) &&
              static_cast<double>(indexIView0(jnode0)-1) <= iMax &&
              jMin <= static_cast<double>(indexJView0(jnode0)-1) &&
              static_cast<double>(indexJView0(jnode0)-1) <= jMax) {
            pointsNeeded = true;
          }
        }
      }
      if (pointsNeeded) {
        ++mRecvCounts_[mpiTask_[i*ny_+j]];
        mRecvPointsList.push_back(i*ny_+j);
      }
    }
  }

  // Buffer size
  mRecvSize_ = mRecvPointsList.size();

  // RecvDispls
  mRecvDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    mRecvDispls_.push_back(mRecvDispls_[jt]+mRecvCounts_[jt]);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> mRecvCountsBuffer(comm_.size());
  comm_.allGatherv(mRecvCounts_.begin(), mRecvCounts_.end(), mRecvCountsBuffer);
  std::vector<int> mRecvCountsGlb_ = std::move(mRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rSendCounts_.push_back(mRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  rSendSize_ = 0;
  for (const auto & n : rSendCounts_) rSendSize_ += n;

  // SendDispls
  rSendDispls_.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    rSendDispls_.push_back(rSendDispls_[jt]+rSendCounts_[jt]);
  }

  // Ordered received points list
  std::vector<size_t> mRecvOffset(comm_.size(), 0);
  std::vector<int> mRecvPointsListOrdered(mRecvSize_);
  for (size_t jr = 0; jr < mRecvSize_; ++jr) {
    size_t jt = mpiTask_[mRecvPointsList[jr]];
    size_t jro = mRecvDispls_[jt]+mRecvOffset[jt];
    mRecvPointsListOrdered[jro] = mRecvPointsList[jr];
    ++mRecvOffset[jt];
  }
  std::vector<int> rSentPointsList(rSendSize_);
  comm_.allToAllv(mRecvPointsListOrdered.data(), mRecvCounts_.data(), mRecvDispls_.data(),
                  rSentPointsList.data(), rSendCounts_.data(), rSendDispls_.data());

  // Mapping for sent points
  rSendMapping_.resize(rSendSize_);
  for (size_t js = 0; js < rSendSize_; ++js) {
    bool found = false;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      if (static_cast<size_t>(rSentPointsList[js]) ==
        (indexIView(jnode)-1)*ny_+indexJView(jnode)-1) {
        ASSERT(!found);
        rSendMapping_[js] = jnode;
        found = true;
      }
    }
    ASSERT(found);
  }

  // Compute horizontal interpolation
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    // Interpolation element default values
    std::string interpType = "n";
    size_t index1 = nx_;
    size_t index2 = ny_;
    std::vector<std::pair<size_t, double>> operations;
    if (ghostView(jnode0) == 0) {
      // Model grid indices
      const size_t i0 = indexIView0(jnode0)-1;
      const size_t j0 = indexJView0(jnode0)-1;
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
        while ((!foundIJ) && jro < mRecvSize_) {
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j) {
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
        while ((!(foundIJ && foundIpJ)) && jro < mRecvSize_) {
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == (i+1)*ny_+j) {
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
        while ((!(foundIJ && foundIJp)) && jro < mRecvSize_) {
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j+1) {
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
        while ((!(foundIJ && foundIpJ && foundIJp && foundIpJp)) && jro < mRecvSize_) {
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j) {
            ASSERT(!foundIJ);
            jroIJ = jro;
            foundIJ = true;
          }
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == (i+1)*ny_+j) {
            ASSERT(!foundIpJ);
            jroIpJ = jro;
            foundIpJ = true;
          }
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == i*ny_+j+1) {
            ASSERT(!foundIJp);
            jroIJp = jro;
            foundIJp = true;
          }
          if (static_cast<size_t>(mRecvPointsListOrdered[jro]) == (i+1)*ny_+j+1) {
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
        if (std::abs(normVertCoord_[k0] - zCoord[k]) < 1.0e-12) {
          // Colocated point
          index1 = k;
          operations.push_back(std::make_pair(k, 1.0));
          found = true;
        }
        ++k;
      }
      k = 0;
      while ((!found) && (k < nz_-1)) {
        if (zCoord[k] < normVertCoord_[k0] && normVertCoord_[k0] < zCoord[k+1]) {
          // Linear interpolation
          index1 = k;
          const double alphaK = normVertCoord_[k0]-zCoord[k];
          operations.push_back(std::make_pair(k, 1.0-alphaK));
          operations.push_back(std::make_pair(k+1, alphaK));
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
    atlas::Field redField = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field modelField = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    auto redView = atlas::array::make_view<double, 2>(redField);
    auto modelView = atlas::array::make_view<double, 2>(modelField);
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      const double x = static_cast<double>(indexIView(jnode)-1)/static_cast<double>(nx_-1);
      const double y = static_cast<double>(indexJView(jnode)-1)/static_cast<double>(ny_-1);
      for (size_t k = 0; k < nz_; ++k) {
        const double z = zCoord[k]/static_cast<double>(nz0_-1);
        redView(jnode, k) = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)
          *std::cos(2.0*M_PI*z)+1.0);
      }
    }
    interpolationTL(redField, modelField);
    double accuracy = 0.0;
    double maxVal = 0.0;
    double maxRefVal = 0.0;
    std::vector<double> locMax(3, 0.0);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        const double x = static_cast<double>(indexIView0(jnode0)-1)/static_cast<double>(nx0_-1);
        const double y = static_cast<double>(indexJView0(jnode0)-1)/static_cast<double>(ny0_-1);
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          const double z = normVertCoord_[k0]/static_cast<double>(nz0_-1);
          double refVal = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)*std::cos(2.0*M_PI*z)+1.0);
          double diff = std::abs(modelView(jnode0, k0)-refVal);
          if (diff > accuracy) {
            accuracy = diff;
            maxVal = modelView(jnode0, k0);
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
    atlas::Field redFieldTL = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field modelFieldTL = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field modelFieldAD = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    atlas::Field redFieldAD = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    auto redViewTL = atlas::array::make_view<double, 2>(redFieldTL);
    auto modelViewTL = atlas::array::make_view<double, 2>(modelFieldTL);
    auto modelViewAD = atlas::array::make_view<double, 2>(modelFieldAD);
    auto redViewAD = atlas::array::make_view<double, 2>(redFieldAD);

    // Generate random fields
    size_t seed = 7;  // To avoid impact on future random generator calls
    util::NormalDistribution<double> dist(rSize_*nz_+mSize_*nz0_, 0.0, 1.0, seed);
    size_t jj = 0;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        redViewTL(jnode, k) = dist[jj];
        ++jj;
      }
    }
    modelViewAD.assign(0.0);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          modelViewAD(jnode0, k) = dist[jj];
          ++jj;
        }
      }
    }

    // Interpolation TL/AD
    interpolationTL(redFieldTL, modelFieldTL);
    interpolationAD(modelFieldAD, redFieldAD);

    // Adjoint test
    double dp1 = 0.0;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      for (size_t k = 0; k < nz_; ++k) {
        dp1 += redViewTL(jnode, k)*redViewAD(jnode, k);
      }
    }
    comm_.allReduceInPlace(dp1, eckit::mpi::sum());
    double dp2 = 0.0;
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          dp2 += modelViewTL(jnode0, k0)*modelViewAD(jnode0, k0);
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

void LayerBase::setupKernels() {
  oops::Log::trace() << classname() << "::setupKernels starting" << std::endl;

  // Get kernels size
  xKernelSize_ = 2*static_cast<size_t>((0.5*rh_+1.0e-12)/xRedFac_)+1;
  yKernelSize_ = 2*static_cast<size_t>((0.5*rh_+1.0e-12)/yRedFac_)+1;
  zKernelSize_ = 2*static_cast<size_t>((0.5*rv_+1.0e-12)/zRedFac_)+1;

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

void LayerBase::read(const int & id) {
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

void LayerBase::broadcast() {
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

std::vector<int> LayerBase::writeDef(const int & id) const {
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

void LayerBase::writeData(const std::vector<int> & varIds) const {
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

void LayerBase::interpolationTL(const atlas::Field & redField,
                                atlas::Field & modelField) const {
  oops::Log::trace() << classname() << "::interpolationTL starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> rSendCounts3D(comm_.size());
  std::vector<int> rSendDispls3D(comm_.size());
  std::vector<int> mRecvCounts3D(comm_.size());
  std::vector<int> mRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rSendCounts3D[jt] = rSendCounts_[jt]*nz_;
    rSendDispls3D[jt] = rSendDispls_[jt]*nz_;
    mRecvCounts3D[jt] = mRecvCounts_[jt]*nz_;
    mRecvDispls3D[jt] = mRecvDispls_[jt]*nz_;
  }

  // Serialize
  const auto redView = atlas::array::make_view<double, 2>(redField);
  std::vector<double> rSendVec(rSendSize_*nz_);
  for (size_t js = 0; js < rSendSize_; ++js) {
    size_t jnode = rSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      ASSERT(js*nz_+k < rSendSize_*nz_);
      rSendVec[js*nz_+k] = redView(jnode, k);
    }
  }

  // Communication
  std::vector<double> mRecvVec(mRecvSize_*nz_);
  comm_.allToAllv(rSendVec.data(), rSendCounts3D.data(), rSendDispls3D.data(),
                  mRecvVec.data(), mRecvCounts3D.data(), mRecvDispls3D.data());

  // Interpolation
  auto modelView = atlas::array::make_view<double, 2>(modelField);
  modelView.assign(0.0);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    for (const auto & horOperation : horInterp_[jnode0].operations()) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        for (const auto & verOperation : verInterp_[k0].operations()) {
          size_t mIndex = horOperation.first*nz_+verOperation.first;
          modelView(jnode0, k0) += horOperation.second*verOperation.second*mRecvVec[mIndex];
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::interpolationTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerBase::interpolationAD(const atlas::Field & modelField,
                                atlas::Field & redField) const {
  oops::Log::trace() << classname() << "::interpolationAD starting" << std::endl;

  // Scale counts and displs for all levels
  std::vector<int> rSendCounts3D(comm_.size());
  std::vector<int> rSendDispls3D(comm_.size());
  std::vector<int> mRecvCounts3D(comm_.size());
  std::vector<int> mRecvDispls3D(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rSendCounts3D[jt] = rSendCounts_[jt]*nz_;
    rSendDispls3D[jt] = rSendDispls_[jt]*nz_;
    mRecvCounts3D[jt] = mRecvCounts_[jt]*nz_;
    mRecvDispls3D[jt] = mRecvDispls_[jt]*nz_;
  }

  // Interpolation adjoint
  const auto modelView = atlas::array::make_view<double, 2>(modelField);
  std::vector<double> mRecvVec(mRecvSize_*nz_, 0.0);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    for (const auto & horOperation : horInterp_[jnode0].operations()) {
      for (size_t k0 = 0; k0 < nz0_; ++k0) {
        for (const auto & verOperation : verInterp_[k0].operations()) {
          size_t mIndex = horOperation.first*nz_+verOperation.first;
          mRecvVec[mIndex] += horOperation.second*verOperation.second*modelView(jnode0, k0);
        }
      }
    }
  }

  // Communication
  std::vector<double> rSendVec(rSendSize_*nz_);
  comm_.allToAllv(mRecvVec.data(), mRecvCounts3D.data(), mRecvDispls3D.data(),
                  rSendVec.data(), rSendCounts3D.data(), rSendDispls3D.data());

  // Deserialize
  auto redView = atlas::array::make_view<double, 2>(redField);
  redView.assign(0.0);
  for (size_t js = 0; js < rSendSize_; ++js) {
    size_t jnode = rSendMapping_[js];
    for (size_t k = 0; k < nz_; ++k) {
      redView(jnode, k) += rSendVec[js*nz_+k];
    }
  }

  oops::Log::trace() << classname() << "::interpolationAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
