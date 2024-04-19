/*
 * (C) Copyright 2024 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/LayerBase.h"

#include <netcdf.h>

#include <algorithm>
#include <unordered_map>
#include <utility>

#include "atlas/array.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "oops/generic/gc99.h"
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
  const eckit::LocalConfiguration & fieldsMetaData,
  const oops::GeometryData & gdata,
  const std::string & myGroup,
  const std::vector<std::string> & myVars,
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
    jsb->second->make(params, fieldsMetaData, gdata, myGroup, myVars, nx0, ny0, nz0);
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
    const std::string key = myGroup_ + ".vert_coord";
    const std::string vertCoordName = fieldsMetaData_.getString(key, "vert_coord");
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          double VC = static_cast<double>(k0+1);
          if (gdata_.fieldSet().has(vertCoordName)) {
            const atlas::Field vertCoordField = gdata_.fieldSet()[vertCoordName];
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

    // Check if vertical length-scale is positive
    bool posRv = true;
    for (size_t k0 = 0; k0 < nz0_; ++k0) {
      if (rv[k0] == 0.0) {
        posRv = false;
        break;
      }
    }

    if (posRv) {
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
    } else {
      normVertCoord_.resize(nz0_, 0.0);
      rv_ = 0.0;
    }
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

  // Interpolation flag
  noInterp_ = !(xRedFac_ > 1.0 || yRedFac_ > 1.0 || zRedFac_ > 1.0);

  // Reduced grid coordinates
  std::vector<double> xCoord;
  const double dxCoord = static_cast<double>(nx0_-1)/static_cast<double>(nx_-1);
  for (size_t i = 0; i < nx_; ++i) {
    xCoord.push_back(static_cast<double>(i)*dxCoord);
  }
  std::vector<double> yCoord;
  const double dyCoord = static_cast<double>(ny0_-1)/static_cast<double>(ny_-1);
  for (size_t j = 0; j < ny_; ++j) {
    yCoord.push_back(static_cast<double>(j)*dyCoord);
  }
  std::vector<double> zCoord;
  if (nz_ > 1) {
    const double dz = static_cast<double>(nz0_-1)/static_cast<double>(nz_-1);
    for (size_t k = 0; k < nz_; ++k) {
      zCoord.push_back(static_cast<double>(k)*dz);
    }
  } else {
    zCoord.push_back(0.0);
  }

  if (noInterp_) {
    if (params_.parallelization.value() == "halo") {
      // Define reduced grid horizontal distribution
      mpiTask_.resize(nx_*ny_, 0);
      std::vector<int> mpiMask(nx_*ny_, 0);
      for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
        if (ghostView(jnode0) == 0) {
          int srcI = indexIView0(jnode0)-1;
          int srcJ = indexJView0(jnode0)-1;
          mpiTask_[srcI*ny_+srcJ] = myrank_;
          mpiMask[srcI*ny_+srcJ] = 1;
        }
      }

      // Check that every point is assigned to a task
      comm_.allReduceInPlace(mpiTask_.begin(), mpiTask_.end(), eckit::mpi::sum());
      comm_.allReduceInPlace(mpiMask.begin(), mpiMask.end(), eckit::mpi::sum());
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t i = 0; i < nx_; ++i) {
          if (mpiMask[i*ny_+j] == 0) {
            oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")"
              << std::endl;
            throw eckit::Exception("task not define for this point", Here());
          }
          if (mpiMask[i*ny_+j] > 1) {
            oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")"
              << std::endl;
            throw eckit::Exception("task defined more than once for this point", Here());
          }
        }
      }
    }

    // Create reduced grid FunctionSpace on each task
    std::vector<atlas::PointXY> v;
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        // Fake coordinate in [0,1]
        const size_t i = indexIView0(jnode0)-1;
        const size_t j = indexJView0(jnode0)-1;
        atlas::PointXY p({xCoord[i]/static_cast<double>(nx0_-1),
          yCoord[j]/static_cast<double>(ny0_-1)});
        v.push_back(p);
      }
    }
    rSize_ = v.size();
    fspace_ = atlas::functionspace::PointCloud(v);

    // Create reduced grid index fields
    atlas::Field fieldIndexI = fspace_.createField<int>(atlas::option::name("index_i"));
    atlas::Field fieldIndexJ = fspace_.createField<int>(atlas::option::name("index_j"));
    auto indexIView = atlas::array::make_view<int, 1>(fieldIndexI);
    auto indexJView = atlas::array::make_view<int, 1>(fieldIndexJ);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        indexIView(jnode0) = indexIView0(jnode0);
        indexJView(jnode0) = indexJView0(jnode0);
      }
    }
    fset_.add(fieldIndexI);
    fset_.add(fieldIndexJ);
  } else {
    // Define reduced grid horizontal distribution
    mpiTask_.resize(nx_*ny_, 0);
    std::vector<int> mpiMask(nx_*ny_, 0);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        int srcI = -1;
        for (size_t i = 0; i < nx_; ++i) {
          if (indexIView0(jnode0)-1 == std::round(xCoord[i])) {
            srcI = i;
            break;
          }
        }
        int srcJ = -1;
        for (size_t j = 0; j < ny_; ++j) {
          if (indexJView0(jnode0)-1 == std::round(yCoord[j])) {
            srcJ = j;
            break;
          }
        }
        if ((srcI > -1) && (srcJ > -1)) {
          mpiTask_[srcI*ny_+srcJ] = myrank_;
          mpiMask[srcI*ny_+srcJ] = 1;
        }
      }
    }

    // Check that every point is assigned to a task
    comm_.allReduceInPlace(mpiTask_.begin(), mpiTask_.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(mpiMask.begin(), mpiMask.end(), eckit::mpi::sum());
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t i = 0; i < nx_; ++i) {
        if (mpiMask[i*ny_+j] == 0) {
          oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")"
            << std::endl;
          throw eckit::Exception("task not define for this point", Here());
        }
        if (mpiMask[i*ny_+j] > 1) {
          oops::Log::info() << "Info     :     Point (i,j) = (" << i << "," << j << ")"
            << std::endl;
          throw eckit::Exception("task defined more than once for this point", Here());
        }
      }
    }

    // Create reduced grid FunctionSpace on each task
    std::vector<atlas::PointXY> v;
    for (size_t j = 0; j < ny_; ++j) {
      for (size_t i = 0; i < nx_; ++i) {
        if (static_cast<size_t>(mpiTask_[i*ny_+j]) == myrank_) {
          // Fake coordinate in [0,1]
          atlas::PointXY p({xCoord[i]/static_cast<double>(nx0_-1),
            yCoord[j]/static_cast<double>(ny0_-1)});
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

    // Define local tree on model grid
    std::vector<atlas::Point3> points0(mSize_);
    std::vector<size_t> indices0(mSize_);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      const double x0 = static_cast<double>(indexIView0(jnode0)-1);
      const double y0 = static_cast<double>(indexJView0(jnode0)-1);
      points0[jnode0] = atlas::Point3(x0, y0, 0.0);
      indices0[jnode0] = jnode0;
    }
    atlas::util::IndexKDTree mTree;
    mTree.build(points0, indices0);
    const double radius = std::sqrt(dxCoord*dxCoord+dyCoord*dyCoord);

    // RecvCounts and received points list
    mRecvCounts_.resize(comm_.size());
    std::fill(mRecvCounts_.begin(), mRecvCounts_.end(), 0);
    std::vector<int> mRecvPointsList;
    for (size_t j = 0; j < ny_; ++j) {
      const double jMin = j > 0 ? yCoord[j-1] : yCoord[0];
      const double jMax = yCoord[std::min(j+1, ny_-1)];
      for (size_t i = 0; i < nx_; ++i) {
        const double iMin = i > 0 ? xCoord[i-1] : xCoord[0];
        const double iMax = xCoord[std::min(i+1, nx_-1)];
        const atlas::Point3 p(xCoord[i], yCoord[j], 0.0);
        const auto list = mTree.closestPointsWithinRadius(p, radius);
        bool pointsNeeded = false;
        for (const auto & item : list) {
          const size_t jnode0 = item.payload();
          if (ghostView(jnode0) == 0) {
            if (iMin <= static_cast<double>(indexIView0(jnode0)-1) &&
              static_cast<double>(indexIView0(jnode0)-1) <= iMax &&
              jMin <= static_cast<double>(indexJView0(jnode0)-1) &&
              static_cast<double>(indexJView0(jnode0)-1) <= jMax) {
              pointsNeeded = true;
              break;
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

    // Sort indices
    std::vector<size_t> gij;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      gij.push_back((indexIView(jnode)-1)*ny_+indexJView(jnode)-1);
    }
    std::vector<size_t> gidx(rSize_);
    std::iota(gidx.begin(), gidx.end(), 0);
    std::stable_sort(gidx.begin(), gidx.end(),
      [&gij](size_t i1, size_t i2) {return gij[i1] < gij[i2];});
    std::vector<size_t> ridx(rSendSize_);
    std::iota(ridx.begin(), ridx.end(), 0);
    std::stable_sort(ridx.begin(), ridx.end(),
      [&rSentPointsList](size_t i1, size_t i2) {return rSentPointsList[i1] < rSentPointsList[i2];});

    // Mapping for sent points
    rSendMapping_.resize(rSendSize_);
    jnode = 0;
    for (size_t js = 0; js < rSendSize_; ++js) {
      while (gij[gidx[jnode]] < static_cast<size_t>(rSentPointsList[ridx[js]])) {
        ++jnode;
        ASSERT(jnode < rSize_);
      }
      rSendMapping_[ridx[js]] = gidx[jnode];
    }

    // Sort indices
    std::vector<size_t> idx(mRecvPointsListOrdered.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&mRecvPointsListOrdered](size_t i1, size_t i2)
      {return mRecvPointsListOrdered[i1] < mRecvPointsListOrdered[i2];});

    // Compute horizontal interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        // Interpolation element default values
        std::string interpType = "n";
        size_t indexI = nx_;
        size_t indexJ = ny_;
        std::vector<std::pair<size_t, double>> operations;

        // Model grid indices
        const size_t i0 = indexIView0(jnode0)-1;
        const size_t j0 = indexJView0(jnode0)-1;
        const double di = static_cast<double>(i0)/xRedFac_;
        const double dj = static_cast<double>(j0)/yRedFac_;
        indexI = static_cast<size_t>(di);
        indexJ = static_cast<size_t>(dj);
        const bool colocatedI = (std::abs(static_cast<double>(indexI)-di) < 1.0e-8);
        const bool colocatedJ = (std::abs(static_cast<double>(indexJ)-dj) < 1.0e-8);
        const double alphaI = di-static_cast<double>(indexI);
        const double alphaJ = dj-static_cast<double>(indexJ);

        // Points to find
        std::vector<bool> toFind = {true, !colocatedI, !colocatedJ, !colocatedI && !colocatedJ};
        std::vector<size_t> valueToFind = {indexI*ny_+indexJ, (indexI+1)*ny_+indexJ,
          indexI*ny_+(indexJ+1), (indexI+1)*ny_+(indexJ+1)};
        std::vector<int> foundIndex(4, -1);

        // Binary search for each point
        for (size_t jj = 0; jj < 4; ++jj) {
          if (toFind[jj]) {
            binarySearch(mRecvPointsListOrdered, idx, valueToFind[jj], foundIndex[jj]);
            ASSERT(foundIndex[jj] > -1);
            ASSERT(static_cast<size_t>(mRecvPointsListOrdered[foundIndex[jj]]) == valueToFind[jj]);
          }
        }

        // Create interpolation operations
        if (colocatedI && colocatedJ) {
          // Colocated point
          interpType = "c";
          operations.push_back(std::make_pair(foundIndex[0], 1.0));
        } else if (colocatedJ) {
          // Linear interpolation along x
          interpType = "x";
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaI));
          operations.push_back(std::make_pair(foundIndex[1], alphaI));
        } else if (colocatedI) {
          // Linear interpolation along y
          interpType = "y";
          operations.push_back(std::make_pair(foundIndex[0], 1.0-alphaJ));
          operations.push_back(std::make_pair(foundIndex[2], alphaJ));
        } else {
          // Bilinear interpolation
          interpType = "b";
          operations.push_back(std::make_pair(foundIndex[0], (1.0-alphaI)*(1.0-alphaJ)));
          operations.push_back(std::make_pair(foundIndex[1], alphaI*(1.0-alphaJ)));
          operations.push_back(std::make_pair(foundIndex[2], (1.0-alphaI)*alphaJ));
          operations.push_back(std::make_pair(foundIndex[3], alphaI*alphaJ));
        }
        horInterp_.push_back(InterpElement(interpType, indexI, indexJ, operations));
      }
    }
  }

  // Compute vertical interpolation
  for (size_t k0 = 0; k0 < nz0_; ++k0) {
    size_t index1 = nz_;
    std::vector<std::pair<size_t, double>> operations;
    if (rv_ > 0.0) {
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
        size_t k = 1;
        while ((!found) && (k < nz_-1)) {
          if (std::abs(normVertCoord_[k0]-zCoord[k]) < 1.0e-12) {
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
    } else {
      // No interpolation
      index1 = k0;
      operations.push_back(std::make_pair(k0, 1.0));
    }
    verInterp_.push_back(InterpElement(index1, operations));
  }

  if (!params_.skipTests.value()) {
    // Get indices
    const auto indexIView = atlas::array::make_view<int, 1>(fset_["index_i"]);
    const auto indexJView = atlas::array::make_view<int, 1>(fset_["index_j"]);

    // Check interpolation accuracy
    atlas::Field redField = fspace_.createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz_));
    atlas::Field modelField = gdata_.functionSpace().createField<double>(
      atlas::option::name("dummy") | atlas::option::levels(nz0_));
    auto redView = atlas::array::make_view<double, 2>(redField);
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
    const auto modelView = atlas::array::make_view<double, 2>(modelField);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        const double x = static_cast<double>(indexIView0(jnode0)-1)/static_cast<double>(nx0_-1);
        const double y = static_cast<double>(indexJView0(jnode0)-1)/static_cast<double>(ny0_-1);
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          const double z = normVertCoord_[k0]/static_cast<double>(nz0_-1);
          const double refVal = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)
            *std::cos(2.0*M_PI*z)+1.0);
          const double diff = std::abs(modelView(jnode0, k0)-refVal);
          if (diff > accuracy) {
            accuracy = diff;
            maxVal = modelView(jnode0, k0);
            maxRefVal = refVal;
            const auto lonLatView = atlas::array::make_view<double, 2>(
              gdata_.functionSpace().lonlat());
            locMax = {lonLatView(jnode0, 0), lonLatView(jnode0, 1), z};
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
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          modelViewAD(jnode0, k0) = dist[jj];
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
  double zNorm = 0.0;
  if (rv_ > 0.0) {
    double zAlpha = zRedFac_/(0.5*rv_);
    for (size_t jk = 0; jk < zKernelSize_; ++jk) {
      int jkc = jk-(zKernelSize_-1)/2;
      if (jkc < 0) {
        zKernel_[jk] = zAlpha*static_cast<double>(jkc)+1.0;
      } else {
        zKernel_[jk] = -zAlpha*static_cast<double>(jkc)+1.0;
      }
      zNorm += zKernel_[jk]*zKernel_[jk];
    }
  } else {
    zKernel_[0] = 1.0;
    zNorm = 1.0;
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
  if (zNorm > 0.0) {
    zNorm = 1.0/std::sqrt(zNorm);
    for (size_t jk = 0; jk < zKernelSize_; ++jk) {
      zKernel_[jk] *= zNorm;
    }
  }

  oops::Log::trace() << classname() << "::setupKernels done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerBase::setupNormalization() {
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

  // Cost-efficient normalization

  // Extract convolution values
  size_t nxHalf = xNormSize_+1;
  size_t nyHalf = yNormSize_+1;
  std::vector<double> horConv(4*nxHalf*nyHalf, 0.0);
  std::vector<double> verConv(2*(nz_-1), 0.0);
  extractConvolution(nxHalf, nyHalf, horConv, verConv);

  // Full normalization
  atlas::Field normField = gdata_.functionSpace().createField<double>(
    atlas::option::name(myGroup_) | atlas::option::levels(nz0_));
  auto normView = atlas::array::make_view<double, 2>(normField);
  norm_.add(normField);

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Compute horizontal normalization
  normView.assign(0.0);
  if (noInterp_) {
    // No horizontal interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        normView(jnode0, 0) = 1.0;
      }
    }
  } else {
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        // Define offset
        size_t offsetI = std::min(std::min(horInterp_[jnode0].index1(),
          nx_-2-horInterp_[jnode0].index1()), nxHalf-1);
        size_t offsetJ = std::min(std::min(horInterp_[jnode0].index2(),
          ny_-2-horInterp_[jnode0].index2()), nyHalf-1);
        size_t horOffset = 4*(offsetI*nyHalf+offsetJ);

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
  }

  // Compute vertical normalization
  std::vector<double> verNorm(nz0_, 1.0);
  if (nz0_ > 1) {
    verNorm[0] = 1.0;
    verNorm[nz0_-1] = 1.0;
    for (size_t k0 = 1; k0 < nz0_-1; ++k0) {
      if (verInterp_[k0].operations().size() > 1) {
        size_t verOffset = 2*verInterp_[k0].index1();
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
    atlas::Field cv("genericCtlVec", atlas::array::make_datatype<double>(),
      atlas::array::make_shape(ctlVecSize()));
    auto modelView = atlas::array::make_view<double, 2>(modelField);

    // Model grid indices
    atlas::Field fieldIndexI0 = gdata_.fieldSet()["index_i"];
    atlas::Field fieldIndexJ0 = gdata_.fieldSet()["index_j"];
    auto indexIView0 = atlas::array::make_view<int, 1>(fieldIndexI0);
    auto indexJView0 = atlas::array::make_view<int, 1>(fieldIndexJ0);

    // Compute normalization accuracy
    oops::Log::info() << "Info     :     Compute exact normalization" << std::endl;
    atlas::Field normAccField = gdata_.functionSpace().createField<double>(
      atlas::option::name(myGroup_) | atlas::option::levels(nz0_));
    auto normAccView = atlas::array::make_view<double, 2>(normAccField);
    normAccView.assign(util::missingValue<double>());
    double normAccMax = 0.0;;

    // Sort indices
    std::vector<int> gij0;
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      gij0.push_back((indexIView0(jnode0)-1)*ny0_+indexJView0(jnode0)-1);
    }
    std::vector<size_t> gidx(mSize_);
    std::iota(gidx.begin(), gidx.end(), 0);
    std::stable_sort(gidx.begin(), gidx.end(),
      [&gij0](size_t i1, size_t i2) {return gij0[i1] < gij0[i2];});

    for (size_t i0 = 0; i0 < nx0_; i0 += params_.normAccStride.value()) {
      for (size_t j0 = 0; j0 < ny0_; j0 += params_.normAccStride.value()) {
        // Binary search
        size_t valueToFind = i0*ny0_+j0;
        int myJnode0;
        binarySearch(gij0, gidx, valueToFind, myJnode0);

        for (size_t k0 = 0; k0 < nz0_; k0 += params_.normAccStride.value()) {
          // Set Dirac point
          modelView.assign(0.0);
          if (myJnode0 > -1) {
            modelView(myJnode0, k0) = 1.0;
          }

          // Adjoint square-root multiplication
          const size_t offset = 0;
          multiplySqrtTrans(modelField, cv, offset);

          // Compute exact normalization
          double exactNorm = 0.0;
          const auto cvView = atlas::array::make_view<double, 1>(cv);
          for (size_t jj = 0; jj < ctlVecSize(); ++jj) {
            exactNorm += cvView(jj)*cvView(jj);
          }
          comm_.allReduceInPlace(exactNorm, eckit::mpi::sum());

          // Get exact normalization factor
          ASSERT(exactNorm > 0.0);
          exactNorm = 1.0/std::sqrt(exactNorm);

          // Assess cost-efficient normalization quality
          if (myJnode0 > -1) {
            normAccView(myJnode0, k0) = (normView(myJnode0, k0)-exactNorm)/exactNorm;
            normAccMax = std::max(normAccMax, std::abs(normAccView(myJnode0, k0)));
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
    if (*params_.resol.value() != resol_) {
      throw eckit::UserError("resol parameter inconsistent between yaml and file", Here());
    }
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

  // Initialization
  const auto redView = atlas::array::make_view<double, 2>(redField);
  auto modelView = atlas::array::make_view<double, 2>(modelField);
  modelView.assign(0.0);

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  if (noInterp_) {
    // No interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          for (const auto & verOperation : verInterp_[k0].operations()) {
            modelView(jnode0, k0) += verOperation.second*redView(jnode0, verOperation.first);
          }
        }
      }
    }
  } else {
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
    std::vector<double> rSendVec(rSendSize_*nz_);
    for (size_t js = 0; js < rSendSize_; ++js) {
      const size_t jnode = rSendMapping_[js];
      for (size_t k = 0; k < nz_; ++k) {
        rSendVec[js*nz_+k] = redView(jnode, k);
      }
    }

    // Communication
    std::vector<double> mRecvVec(mRecvSize_*nz_);
    comm_.allToAllv(rSendVec.data(), rSendCounts3D.data(), rSendDispls3D.data(),
      mRecvVec.data(), mRecvCounts3D.data(), mRecvDispls3D.data());

    // Interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (const auto & horOperation : horInterp_[jnode0].operations()) {
          for (size_t k0 = 0; k0 < nz0_; ++k0) {
            for (const auto & verOperation : verInterp_[k0].operations()) {
              const size_t mIndex = horOperation.first*nz_+verOperation.first;
              modelView(jnode0, k0) += horOperation.second*verOperation.second*mRecvVec[mIndex];
            }
          }
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

  // Initialization
  const auto modelView = atlas::array::make_view<double, 2>(modelField);
  auto redView = atlas::array::make_view<double, 2>(redField);
  redView.assign(0.0);

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  if (noInterp_) {
    // No interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t k0 = 0; k0 < nz0_; ++k0) {
          for (const auto & verOperation : verInterp_[k0].operations()) {
            redView(jnode0, verOperation.first) += verOperation.second*modelView(jnode0, k0);
          }
        }
      }
    }
  } else {
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
    std::vector<double> mRecvVec(mRecvSize_*nz_, 0.0);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (const auto & horOperation : horInterp_[jnode0].operations()) {
          for (size_t k0 = 0; k0 < nz0_; ++k0) {
            for (const auto & verOperation : verInterp_[k0].operations()) {
              const size_t mIndex = horOperation.first*nz_+verOperation.first;
              mRecvVec[mIndex] += horOperation.second*verOperation.second*modelView(jnode0, k0);
            }
          }
        }
      }
    }

    // Communication
    std::vector<double> rSendVec(rSendSize_*nz_);
    comm_.allToAllv(mRecvVec.data(), mRecvCounts3D.data(), mRecvDispls3D.data(),
      rSendVec.data(), rSendCounts3D.data(), rSendDispls3D.data());

    // Deserialize
    for (size_t js = 0; js < rSendSize_; ++js) {
      const size_t jnode = rSendMapping_[js];
      for (size_t k = 0; k < nz_; ++k) {
        redView(jnode, k) += rSendVec[js*nz_+k];
      }
    }
  }

  oops::Log::trace() << classname() << "::interpolationAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerBase::binarySearch(const std::vector<int> & vec,
                             const std::vector<size_t> & idx,
                             const size_t & valueToFind,
                             int & foundIndex) {
  // Initialization
  foundIndex = -1;
  size_t low = 0;
  size_t high = vec.size()-1;

  // Loop
  while (low <= high) {
    size_t mid = low+(high-low)/2;
    if (valueToFind == static_cast<size_t>(vec[idx[mid]])) {
      foundIndex = idx[mid];
      break;
    }
    if (valueToFind > static_cast<size_t>(vec[idx[mid]])) {
      low = mid+1;
    }
    if (valueToFind < static_cast<size_t>(vec[idx[mid]])) {
      if (mid == 0) {
        break;
      }
      high = mid-1;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
