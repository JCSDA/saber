/*
 * (C) Copyright 2024 Meteorologisk Institutt
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
        for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
          double VC = static_cast<double>(jz0+1);
          if (gdata_.fieldSet().has(vertCoordName)) {
            const atlas::Field vertCoordField = gdata_.fieldSet()[vertCoordName];
            const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoordField);
            VC = vertCoordView(jnode0, jz0);
          }
          vertCoord[jz0] += VC*wgtView(jnode0, jz0);
          rv[jz0] += rvView(jnode0, jz0)*wgtView(jnode0, jz0);
          wgt[jz0] += wgtView(jnode0, jz0);
        }
      }
    }
    comm_.allReduceInPlace(vertCoord.begin(), vertCoord.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(rv.begin(), rv.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(wgt.begin(), wgt.end(), eckit::mpi::sum());

    // Apply weight
    for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
      ASSERT(wgt[jz0] > 0.0);
      vertCoord[jz0] = vertCoord[jz0]/wgt[jz0];
      rv[jz0] = rv[jz0]/wgt[jz0];
    }

    // Check if vertical length-scale is positive
    bool posRv = true;
    for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
      if (rv[jz0] == 0.0) {
        posRv = false;
        break;
      }
    }

    if (posRv) {
      // Compute thickness
      std::vector<double> thickness(nz0_, 0.0);
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        if (jz0 == 0) {
          thickness[jz0] = std::abs(vertCoord[jz0+1]-vertCoord[jz0]);
        } else if (jz0 == nz0_-1) {
          thickness[jz0] = std::abs(vertCoord[jz0]-vertCoord[jz0-1]);
        } else {
          thickness[jz0] = 0.5*std::abs(vertCoord[jz0+1]-vertCoord[jz0-1]);
        }
      }

      // Normalize thickness with vertical length-scale
      std::vector<double> normThickness(nz0_, 0.0);
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        ASSERT(rv[jz0] > 0.0);
        normThickness[jz0] = thickness[jz0]/rv[jz0];
      }

      // Compute normalized vertical coordinate
      normVertCoord_.resize(nz0_, 0.0);
      for (size_t jz0 = 1; jz0 < nz0_; ++jz0) {
        normVertCoord_[jz0] = normVertCoord_[jz0-1]+0.5*(normThickness[jz0]+normThickness[jz0-1]);
      }

      // Rescale normalized vertical coordinate from 0 to nz0_-1
      const double maxNormVertCoord = normVertCoord_[nz0_-1];
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        normVertCoord_[jz0] = normVertCoord_[jz0]/maxNormVertCoord*static_cast<double>(nz0_-1);
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
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  atlas::Field indexX0Field = fs.index_i();
  atlas::Field indexY0Field = fs.index_j();
  auto indexX0View = atlas::array::make_view<int, 1>(indexX0Field);
  auto indexY0View = atlas::array::make_view<int, 1>(indexY0Field);

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
  std::vector<double> xCoord(nx_);
  const double dxCoord = static_cast<double>(nx0_-1)/static_cast<double>(nx_-1);
  for (size_t jx = 0; jx < nx_; ++jx) {
    xCoord[jx] = static_cast<double>(jx)*dxCoord;
  }
  std::vector<double> yCoord(ny_);
  const double dyCoord = static_cast<double>(ny0_-1)/static_cast<double>(ny_-1);
  for (size_t jy = 0; jy < ny_; ++jy) {
    yCoord[jy] = static_cast<double>(jy)*dyCoord;
  }
  std::vector<double> zCoord;
  if (nz_ > 1) {
    zCoord.resize(nz_);
    const double dz = static_cast<double>(nz0_-1)/static_cast<double>(nz_-1);
    for (size_t jz = 0; jz < nz_; ++jz) {
      zCoord[jz] = static_cast<double>(jz)*dz;
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
          const int srcX = indexX0View(jnode0)-1;
          const int srcY = indexY0View(jnode0)-1;
          mpiTask_[srcX*ny_+srcY] = myrank_;
          mpiMask[srcX*ny_+srcY] = 1;
        }
      }

      // Check that every point is assigned to a task
      comm_.allReduceInPlace(mpiTask_.begin(), mpiTask_.end(), eckit::mpi::sum());
      comm_.allReduceInPlace(mpiMask.begin(), mpiMask.end(), eckit::mpi::sum());
      for (size_t jy = 0; jy < ny_; ++jy) {
        for (size_t jx = 0; jx < nx_; ++jx) {
          if (mpiMask[jx*ny_+jy] == 0) {
            oops::Log::info() << "Info     :     Point (jx,jy) = (" << jx << "," << jy << ")"
              << std::endl;
            throw eckit::Exception("task not define for this point", Here());
          }
          if (mpiMask[jx*ny_+jy] > 1) {
            oops::Log::info() << "Info     :     Point (jx,jy) = (" << jx << "," << jy << ")"
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
        const size_t jx = indexX0View(jnode0)-1;
        const size_t jy = indexY0View(jnode0)-1;
        atlas::PointXY p({xCoord[jx]/static_cast<double>(nx0_-1),
          yCoord[jy]/static_cast<double>(ny0_-1)});
        v.push_back(p);
      }
    }
    rSize_ = v.size();
    fspace_ = atlas::functionspace::PointCloud(v);

    // Create reduced grid index fields
    atlas::Field indexXField = fspace_.createField<int>(atlas::option::name("indexX"));
    atlas::Field indexYField = fspace_.createField<int>(atlas::option::name("indexY"));
    auto indexXView = atlas::array::make_view<int, 1>(indexXField);
    auto indexYView = atlas::array::make_view<int, 1>(indexYField);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        indexXView(jnode0) = indexX0View(jnode0);
        indexYView(jnode0) = indexY0View(jnode0);
      }
    }
    fset_.add(indexXField);
    fset_.add(indexYField);
  } else {
    // Define reduced grid horizontal distribution
    mpiTask_.resize(nx_*ny_, 0);
    std::vector<int> mpiMask(nx_*ny_, 0);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        int srcX = -1;
        for (size_t jx = 0; jx < nx_; ++jx) {
          if (indexX0View(jnode0)-1 == std::round(xCoord[jx])) {
            srcX = jx;
            break;
          }
        }
        int srcY = -1;
        for (size_t jy = 0; jy < ny_; ++jy) {
          if (indexY0View(jnode0)-1 == std::round(yCoord[jy])) {
            srcY = jy;
            break;
          }
        }
        if ((srcX > -1) && (srcY > -1)) {
          mpiTask_[srcX*ny_+srcY] = myrank_;
          mpiMask[srcX*ny_+srcY] = 1;
        }
      }
    }

    // Check that every point is assigned to a task
    comm_.allReduceInPlace(mpiTask_.begin(), mpiTask_.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(mpiMask.begin(), mpiMask.end(), eckit::mpi::sum());
    for (size_t jy = 0; jy < ny_; ++jy) {
      for (size_t jx = 0; jx < nx_; ++jx) {
        if (mpiMask[jx*ny_+jy] == 0) {
          oops::Log::info() << "Info     :     Point (jx,jy) = (" << jx << "," << jy << ")"
            << std::endl;
          throw eckit::Exception("task not define for this point", Here());
        }
        if (mpiMask[jx*ny_+jy] > 1) {
          oops::Log::info() << "Info     :     Point (jx,jy) = (" << jx << "," << jy << ")"
            << std::endl;
          throw eckit::Exception("task defined more than once for this point", Here());
        }
      }
    }

    // Create reduced grid FunctionSpace on each task
    std::vector<atlas::PointXY> v;
    for (size_t jy = 0; jy < ny_; ++jy) {
      for (size_t jx = 0; jx < nx_; ++jx) {
        if (static_cast<size_t>(mpiTask_[jx*ny_+jy]) == myrank_) {
          // Fake coordinate in [0,1]
          atlas::PointXY p({xCoord[jx]/static_cast<double>(nx0_-1),
            yCoord[jy]/static_cast<double>(ny0_-1)});
          v.push_back(p);
        }
      }
    }
    rSize_ = v.size();
    fspace_ = atlas::functionspace::PointCloud(v);

    // Create reduced grid index fields
    atlas::Field indexXField = fspace_.createField<int>(atlas::option::name("indexX"));
    atlas::Field indexYField = fspace_.createField<int>(atlas::option::name("indexY"));
    auto indexXView = atlas::array::make_view<int, 1>(indexXField);
    auto indexYView = atlas::array::make_view<int, 1>(indexYField);
    size_t jnode = 0;
    for (size_t jy = 0; jy < ny_; ++jy) {
      for (size_t jx = 0; jx < nx_; ++jx) {
        if (static_cast<size_t>(mpiTask_[jx*ny_+jy]) == myrank_) {
          indexXView(jnode) = jx+1;
          indexYView(jnode) = jy+1;
          ++jnode;
        }
      }
    }
    fset_.add(indexXField);
    fset_.add(indexYField);

    // Define local tree on model grid
    std::vector<atlas::Point3> points0(mSize_);
    std::vector<size_t> indices0(mSize_);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      const double x0 = static_cast<double>(indexX0View(jnode0)-1);
      const double y0 = static_cast<double>(indexY0View(jnode0)-1);
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
    for (size_t jy = 0; jy < ny_; ++jy) {
      const double jMin = jy > 0 ? yCoord[jy-1] : yCoord[0];
      const double jMax = yCoord[std::min(jy+1, ny_-1)];
      for (size_t jx = 0; jx < nx_; ++jx) {
        const double iMin = jx > 0 ? xCoord[jx-1] : xCoord[0];
        const double iMax = xCoord[std::min(jx+1, nx_-1)];
        const atlas::Point3 p(xCoord[jx], yCoord[jy], 0.0);
        const auto list = mTree.closestPointsWithinRadius(p, radius);
        bool pointsNeeded = false;
        for (const auto & item : list) {
          const size_t jnode0 = item.payload();
          if (ghostView(jnode0) == 0) {
            if (iMin <= static_cast<double>(indexX0View(jnode0)-1) &&
              static_cast<double>(indexX0View(jnode0)-1) <= iMax &&
              jMin <= static_cast<double>(indexY0View(jnode0)-1) &&
              static_cast<double>(indexY0View(jnode0)-1) <= jMax) {
              pointsNeeded = true;
              break;
            }
          }
        }
        if (pointsNeeded) {
          ++mRecvCounts_[mpiTask_[jx*ny_+jy]];
          mRecvPointsList.push_back(jx*ny_+jy);
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
      const size_t jt = mpiTask_[mRecvPointsList[jr]];
      const size_t jro = mRecvDispls_[jt]+mRecvOffset[jt];
      mRecvPointsListOrdered[jro] = mRecvPointsList[jr];
      ++mRecvOffset[jt];
    }
    std::vector<int> rSentPointsList(rSendSize_);
    comm_.allToAllv(mRecvPointsListOrdered.data(), mRecvCounts_.data(), mRecvDispls_.data(),
      rSentPointsList.data(), rSendCounts_.data(), rSendDispls_.data());

    // Sort indices
    std::vector<size_t> gij;
    for (size_t jnode = 0; jnode < rSize_; ++jnode) {
      gij.push_back((indexXView(jnode)-1)*ny_+indexYView(jnode)-1);
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
    horType_.resize(mSize_);
    horStencil_.resize(mSize_);
    horWeights_.resize(mSize_);
    horStencilSize_.resize(mSize_);
    horIndexX_.resize(mSize_);
    horIndexY_.resize(mSize_);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        // Model grid indices
        const size_t jx0 = indexX0View(jnode0)-1;
        const size_t jy0 = indexY0View(jnode0)-1;
        const double dX = static_cast<double>(jx0)/xRedFac_;
        const double dY = static_cast<double>(jy0)/yRedFac_;
        const size_t indexX = static_cast<size_t>(dX);
        const size_t indexY = static_cast<size_t>(dY);
        const bool colocatedX = (std::abs(static_cast<double>(indexX)-dX) < 1.0e-8);
        const bool colocatedJ = (std::abs(static_cast<double>(indexY)-dY) < 1.0e-8);
        const double alphaX = dX-static_cast<double>(indexX);
        const double alphaY = dY-static_cast<double>(indexY);

        // Points to find
        std::array<bool, 4> toFind = {true, !colocatedX, !colocatedJ, !colocatedX && !colocatedJ};
        std::array<size_t, 4> valueToFind = {indexX*ny_+indexY, (indexX+1)*ny_+indexY,
          indexX*ny_+(indexY+1), (indexX+1)*ny_+(indexY+1)};
        std::array<int, 4> foundIndex;
        foundIndex.fill(-1);

        // Binary search for each point
        for (size_t jj = 0; jj < 4; ++jj) {
          if (toFind[jj]) {
            binarySearch(mRecvPointsListOrdered, idx, valueToFind[jj], foundIndex[jj]);
            ASSERT(foundIndex[jj] > -1);
            ASSERT(static_cast<size_t>(mRecvPointsListOrdered[foundIndex[jj]]) == valueToFind[jj]);
          }
        }

        // Create interpolation operations
        if (colocatedX && colocatedJ) {
          // Colocated point
          horType_[jnode0] = "c";
          horStencil_[jnode0][0] = foundIndex[0];
          horWeights_[jnode0][0] = 1.0;
          horStencilSize_[jnode0] = 1;
        } else if (colocatedJ) {
          // Linear interpolation along x
          horType_[jnode0] = "x";
          horStencil_[jnode0][0] = foundIndex[0];
          horWeights_[jnode0][0] = 1.0-alphaX;
          horStencil_[jnode0][1] = foundIndex[1];
          horWeights_[jnode0][1] = alphaX;
          horStencilSize_[jnode0] = 2;
        } else if (colocatedX) {
          // Linear interpolation along y
          horType_[jnode0] = "y";
          horStencil_[jnode0][0] = foundIndex[0];
          horWeights_[jnode0][0] = 1.0-alphaY;
          horStencil_[jnode0][1] = foundIndex[2];
          horWeights_[jnode0][1] = alphaY;
          horStencilSize_[jnode0] = 2;
        } else {
          // Bilinear interpolation
          horType_[jnode0] = "b";
          horStencil_[jnode0][0] = foundIndex[0];
          horWeights_[jnode0][0] = (1.0-alphaX)*(1.0-alphaY);
          horStencil_[jnode0][1] = foundIndex[1];
          horWeights_[jnode0][1] = alphaX*(1.0-alphaY);
          horStencil_[jnode0][2] = foundIndex[2];
          horWeights_[jnode0][2] = (1.0-alphaX)*alphaY;
          horStencil_[jnode0][3] = foundIndex[3];
          horWeights_[jnode0][3] = alphaX*alphaY;
          horStencilSize_[jnode0] = 4;
        }
        horIndexX_[jnode0] = indexX;
        horIndexY_[jnode0] = indexY;
      } else {
        // No interpolation
        horStencilSize_[jnode0] = 0;
      }
    }
  }

  // Compute vertical interpolation
  verStencil_.resize(nz0_);
  verWeights_.resize(nz0_);
  verStencilSize_.resize(nz0_);
  verIndex_.resize(nz0_);
  for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
    if (rv_ > 0.0) {
      if (jz0 == 0) {
        // First level
        verStencil_[jz0][0] = 0;
        verWeights_[jz0][0] = 1.0;
        verStencilSize_[jz0] = 1;
        verIndex_[jz0] = 0;
      } else if (jz0 == nz0_-1) {
        // Last level
        verStencil_[jz0][0] = nz_-1;
        verWeights_[jz0][0] = 1.0;
        verStencilSize_[jz0] = 1;
        verIndex_[jz0] = nz_-1;
      } else {
        // Other levels
        bool found = false;
        size_t jz = 1;
        while ((!found) && (jz < nz_-1)) {
          if (std::abs(normVertCoord_[jz0]-zCoord[jz]) < 1.0e-12) {
            // Colocated point
            verStencil_[jz0][0] = jz;
            verWeights_[jz0][0] = 1.0;
            verStencilSize_[jz0] = 1;
            verIndex_[jz0] = jz;
            found = true;
          }
          ++jz;
        }
        jz = 0;
        while ((!found) && (jz < nz_-1)) {
          if (zCoord[jz] < normVertCoord_[jz0] && normVertCoord_[jz0] < zCoord[jz+1]) {
            // Linear interpolation
            const double alphaZ = normVertCoord_[jz0]-zCoord[jz];
            verStencil_[jz0][0] = jz;
            verWeights_[jz0][0] = 1.0-alphaZ;
            verStencil_[jz0][1] = jz+1;
            verWeights_[jz0][1] = alphaZ;
            verStencilSize_[jz0] = 2;
            verIndex_[jz0] = jz;
            found = true;
          }
          ++jz;
        }
        ASSERT(found);
      }
    } else {
      // No interpolation
      verStencil_[jz0][0] = jz0;
      verWeights_[jz0][0] = 1.0;
      verStencilSize_[jz0] = 1;
      verIndex_[jz0] = jz0;
    }
  }

  if (!params_.skipTests.value()) {
    // Test interpolation
    testInterpolation(zCoord);
  }

  oops::Log::trace() << classname() << "::setupInterpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void LayerBase::testInterpolation(const std::vector<double> & zCoord) const {
  oops::Log::trace() << classname() << "::testInterpolation starting" << std::endl;

  // Get indices
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  atlas::Field indexX0Field = fs.index_i();
  atlas::Field indexY0Field = fs.index_j();
  const auto indexXView = atlas::array::make_view<int, 1>(fset_["indexX"]);
  const auto indexYView = atlas::array::make_view<int, 1>(fset_["indexY"]);
  auto indexX0View = atlas::array::make_view<int, 1>(indexX0Field);
  auto indexY0View = atlas::array::make_view<int, 1>(indexY0Field);

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Check interpolation accuracy
  atlas::Field redField = fspace_.createField<double>(
    atlas::option::name("dummy") | atlas::option::levels(nz_));
  atlas::Field modelField = gdata_.functionSpace().createField<double>(
    atlas::option::name("dummy") | atlas::option::levels(nz0_));
  auto redView = atlas::array::make_view<double, 2>(redField);
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    const double x = static_cast<double>(indexXView(jnode)-1)/static_cast<double>(nx_-1);
    const double y = static_cast<double>(indexYView(jnode)-1)/static_cast<double>(ny_-1);
    for (size_t jz = 0; jz < nz_; ++jz) {
      const double z = zCoord[jz]/static_cast<double>(nz0_-1);
      redView(jnode, jz) = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)
        *std::cos(2.0*M_PI*z)+1.0);
    }
  }
  interpolationTL(redField, modelField);
  double accuracy = 0.0;
  double maxVal = 0.0;
  double maxRefVal = 0.0;
  std::array<double, 3> locMax;
  locMax.fill(0.0);
  const auto modelView = atlas::array::make_view<double, 2>(modelField);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      const double x = static_cast<double>(indexX0View(jnode0)-1)/static_cast<double>(nx0_-1);
      const double y = static_cast<double>(indexY0View(jnode0)-1)/static_cast<double>(ny0_-1);
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        const double z = normVertCoord_[jz0]/static_cast<double>(nz0_-1);
        const double refVal = 0.5*(std::sin(2.0*M_PI*x)*std::sin(2.0*M_PI*y)
          *std::cos(2.0*M_PI*z)+1.0);
        const double diff = std::abs(modelView(jnode0, jz0)-refVal);
        if (diff > accuracy) {
          accuracy = diff;
          maxVal = modelView(jnode0, jz0);
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
  comm_.broadcast(locMax.begin(), locMax.end(), maxTask);
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
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      const size_t jnz = jnode*nz_ + jz;
      redViewTL(jnode, jz) = dist[jnz];
    }
  }
  modelViewAD.assign(0.0);
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        const size_t jnz0 = jnode0*nz0_ + jz0;
        modelViewAD(jnode0, jz0) = dist[jnz0];
      }
    }
  }

  // Interpolation TL/AD
  interpolationTL(redFieldTL, modelFieldTL);
  interpolationAD(modelFieldAD, redFieldAD);

  // Adjoint test
  double dp1 = 0.0;
  for (size_t jnode = 0; jnode < rSize_; ++jnode) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      dp1 += redViewTL(jnode, jz)*redViewAD(jnode, jz);
    }
  }
  comm_.allReduceInPlace(dp1, eckit::mpi::sum());
  double dp2 = 0.0;
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        dp2 += modelViewTL(jnode0, jz0)*modelViewAD(jnode0, jz0);
      }
    }
  }
  comm_.allReduceInPlace(dp2, eckit::mpi::sum());

  // Print result
  const double diff = std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2));
  oops::Log::info() << std::setprecision(16)
    << "Info     :     FastLAM interpolation adjoint test: y^t (Ax) = "
    << dp1 << ": x^t (A^t y) = " << dp2 << " : diff = " << diff << ", adjoint tolerance = "
    << params_.adjointTolerance.value() << std::endl;
  oops::Log::test() << "    FastLAM interpolation adjoint test";
  if (diff < params_.adjointTolerance.value()) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    throw eckit::Exception("interpolation adjoint test failed for block FastLAM", Here());
  }

  oops::Log::trace() << classname() << "::testInterpolation done" << std::endl;
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
  for (size_t jkx = 0; jkx < xKernelSize_; ++jkx) {
    int jkc = jkx-(xKernelSize_-1)/2;
    if (jkc < 0) {
      xKernel_[jkx] = xAlpha*static_cast<double>(jkc)+1.0;
    } else {
      xKernel_[jkx] = -xAlpha*static_cast<double>(jkc)+1.0;
    }
    xNorm += xKernel_[jkx]*xKernel_[jkx];
  }
  double yAlpha = yRedFac_/(0.5*rh_);
  double yNorm = 0.0;
  for (size_t jky = 0; jky < yKernelSize_; ++jky) {
    int jkc = jky-(yKernelSize_-1)/2;
    if (jkc < 0) {
      yKernel_[jky] = yAlpha*static_cast<double>(jkc)+1.0;
    } else {
      yKernel_[jky] = -yAlpha*static_cast<double>(jkc)+1.0;
    }
    yNorm += yKernel_[jky]*yKernel_[jky];
  }
  double zNorm = 0.0;
  if (rv_ > 0.0) {
    double zAlpha = zRedFac_/(0.5*rv_);
    for (size_t jkz = 0; jkz < zKernelSize_; ++jkz) {
      int jkc = jkz-(zKernelSize_-1)/2;
      if (jkc < 0) {
        zKernel_[jkz] = zAlpha*static_cast<double>(jkc)+1.0;
      } else {
        zKernel_[jkz] = -zAlpha*static_cast<double>(jkc)+1.0;
      }
      zNorm += zKernel_[jkz]*zKernel_[jkz];
    }
  } else {
    zKernel_[0] = 1.0;
    zNorm = 1.0;
  }

  // Normalize kernels
  xNorm = 1.0/std::sqrt(xNorm);
  for (size_t jkx = 0; jkx < xKernelSize_; ++jkx) {
    xKernel_[jkx] *= xNorm;
  }
  yNorm = 1.0/std::sqrt(yNorm);
  for (size_t jky = 0; jky < yKernelSize_; ++jky) {
    yKernel_[jky] *= yNorm;
  }
  if (zNorm > 0.0) {
    zNorm = 1.0/std::sqrt(zNorm);
    for (size_t jkz = 0; jkz < zKernelSize_; ++jkz) {
      zKernel_[jkz] *= zNorm;
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
  for (size_t jnx = 0; jnx < xNormSize_; ++jnx) {
    for (size_t jkx = xNormSize_-jnx; jkx < xKernelSize_; ++jkx) {
      xNorm_[jnx] += xKernel_[jkx]*xKernel_[jkx];
    }
    xNorm_[jnx] = 1.0/std::sqrt(xNorm_[jnx]);
  }
  for (size_t jny = 0; jny < yNormSize_; ++jny) {
    for (size_t jky = yNormSize_-jny; jky < yKernelSize_; ++jky) {
      yNorm_[jny] += yKernel_[jky]*yKernel_[jky];
    }
    yNorm_[jny] = 1.0/std::sqrt(yNorm_[jny]);
  }
  for (size_t jnz = 0; jnz < zNormSize_; ++jnz) {
    for (size_t jkz = zNormSize_-jnz; jkz < zKernelSize_; ++jkz) {
      zNorm_[jnz] += zKernel_[jkz]*zKernel_[jkz];
    }
    zNorm_[jnz] = 1.0/std::sqrt(zNorm_[jnz]);
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
        const size_t offsetX = std::min(std::min(horIndexX_[jnode0], nx_-2-horIndexX_[jnode0]),
          nxHalf-1);
        const size_t offsetY = std::min(std::min(horIndexY_[jnode0], ny_-2-horIndexY_[jnode0]),
          nyHalf-1);
        const size_t horOffset = 4*(offsetX*nyHalf+offsetY);

        if (horType_[jnode0] == "c") {
          // Colocated point, no normalization needed
          normView(jnode0, 0) = 1.0;
        } else if (horType_[jnode0] == "x") {
          // Linear interpolation along x
          double xW = horConv[horOffset+0]*horWeights_[jnode0][0]
            +horConv[horOffset+1]*horWeights_[jnode0][1];
          double xE = horConv[horOffset+1]*horWeights_[jnode0][0]
            +horConv[horOffset+0]*horWeights_[jnode0][1];
          normView(jnode0, 0) = horWeights_[jnode0][0]*xW
            +horWeights_[jnode0][1]*xE;
        } else if (horType_[jnode0] == "y") {
          // Linear interpolation along y
          double xS = horConv[horOffset+0]*horWeights_[jnode0][0]
            +horConv[horOffset+2]*horWeights_[jnode0][1];
          double xN = horConv[horOffset+2]*horWeights_[jnode0][0]
            +horConv[horOffset+0]*horWeights_[jnode0][1];
          normView(jnode0, 0) = horWeights_[jnode0][0]*xS
            +horWeights_[jnode0][1]*xN;
        } else if (horType_[jnode0] == "b") {
          // Bilinear interpolation
          double xSW = horConv[horOffset+0]*horWeights_[jnode0][0]
            +horConv[horOffset+1]*horWeights_[jnode0][1]
            +horConv[horOffset+2]*horWeights_[jnode0][2]
            +horConv[horOffset+3]*horWeights_[jnode0][3];
          double xSE = horConv[horOffset+1]*horWeights_[jnode0][0]
            +horConv[horOffset+0]*horWeights_[jnode0][1]
            +horConv[horOffset+3]*horWeights_[jnode0][2]
            +horConv[horOffset+2]*horWeights_[jnode0][3];
          double xNW = horConv[horOffset+2]*horWeights_[jnode0][0]
            +horConv[horOffset+3]*horWeights_[jnode0][1]
            +horConv[horOffset+0]*horWeights_[jnode0][2]
            +horConv[horOffset+1]*horWeights_[jnode0][3];
          double xNE = horConv[horOffset+3]*horWeights_[jnode0][0]
            +horConv[horOffset+2]*horWeights_[jnode0][1]
            +horConv[horOffset+1]*horWeights_[jnode0][2]
            +horConv[horOffset+0]*horWeights_[jnode0][3];
          normView(jnode0, 0) = horWeights_[jnode0][0]*xSW
            +horWeights_[jnode0][1]*xSE
            +horWeights_[jnode0][2]*xNW
            +horWeights_[jnode0][3]*xNE;
        } else {
          throw eckit::Exception("wrong interpolation type: " + horType_[jnode0],
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
    for (size_t jz0 = 1; jz0 < nz0_-1; ++jz0) {
      if (verStencilSize_[jz0] > 1) {
        const size_t verOffset = 2*verIndex_[jz0];
        const double xB = verConv[verOffset+0]*verWeights_[jz0][0]
          +verConv[verOffset+1]*verWeights_[jz0][1];
        const double xT = verConv[verOffset+1]*verWeights_[jz0][0]
          +verConv[verOffset+0]*verWeights_[jz0][1];
        verNorm[jz0] = verWeights_[jz0][0]*xB+verWeights_[jz0][1]*xT;
      }
    }
  }

  // Compute 3D normalization
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
        normView(jnode0, jz0) = normView(jnode0, 0)*verNorm[jz0];
      }
    }
  }

  // Get normalization factor
  for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
    for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
      if (normView(jnode0, jz0) > 0.0) {
        normView(jnode0, jz0) = 1.0/std::sqrt(normView(jnode0, jz0));
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
    const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
    atlas::Field indexX0Field = fs.index_i();
    atlas::Field indexY0Field = fs.index_j();
    auto indexX0View = atlas::array::make_view<int, 1>(indexX0Field);
    auto indexY0View = atlas::array::make_view<int, 1>(indexY0Field);

    // Compute normalization accuracy
    oops::Log::info() << "Info     :     Compute exact normalization" << std::endl;
    atlas::Field normAccField = gdata_.functionSpace().createField<double>(
      atlas::option::name(myGroup_) | atlas::option::levels(nz0_));
    auto normAccView = atlas::array::make_view<double, 2>(normAccField);
    normAccView.assign(util::missingValue<double>());
    double normAccMax = 0.0;;

    // Sort indices
    std::vector<int> gij0(mSize_);
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      gij0[jnode0] = (indexX0View(jnode0)-1)*ny0_+indexY0View(jnode0)-1;
    }
    std::vector<size_t> gidx(mSize_);
    std::iota(gidx.begin(), gidx.end(), 0);
    std::stable_sort(gidx.begin(), gidx.end(),
      [&gij0](size_t i1, size_t i2) {return gij0[i1] < gij0[i2];});

    for (size_t jx0 = 0; jx0 < nx0_; jx0 += params_.normAccStride.value()) {
      for (size_t jy0 = 0; jy0 < ny0_; jy0 += params_.normAccStride.value()) {
        // Binary search
        const size_t valueToFind = jx0*ny0_+jy0;
        int myJnode0;
        binarySearch(gij0, gidx, valueToFind, myJnode0);

        for (size_t jz0 = 0; jz0 < nz0_; jz0 += params_.normAccStride.value()) {
          // Set Dirac point
          modelView.assign(0.0);
          if (myJnode0 > -1) {
            modelView(myJnode0, jz0) = 1.0;
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
            normAccView(myJnode0, jz0) = (normView(myJnode0, jz0)-exactNorm)/exactNorm;
            normAccMax = std::max(normAccMax, std::abs(normAccView(myJnode0, jz0)));
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

  // Resize vectors
  normVertCoord_.resize(nz0_);
  xKernel_.resize(xKernelSize_);
  yKernel_.resize(yKernelSize_);
  zKernel_.resize(zKernelSize_);
  xNorm_.resize(xNormSize_);
  yNorm_.resize(yNormSize_);
  zNorm_.resize(zNormSize_);

  // Read data
  if ((retval = nc_get_var_double(id, normVertCoord_id, normVertCoord_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, xKernel_id, xKernel_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, yKernel_id, yKernel_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, zKernel_id, zKernel_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, xNorm_id, xNorm_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, yNorm_id, yNorm_.data()))) ERR(retval);
  if ((retval = nc_get_var_double(id, zNorm_id, zNorm_.data()))) ERR(retval);

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

std::array<int, 8> LayerBase::writeDef(const int & id) const {
  oops::Log::trace() << classname() << "::writeDef starting" << std::endl;

  // Vector of ids
  std::array<int, 8> varIds;

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

void LayerBase::writeData(const std::array<int, 8> & varIds) const {
  oops::Log::trace() << classname() << "::writeData starting" << std::endl;

  // Write data
  int retval;
  if ((retval = nc_put_var_double(varIds[0], varIds[1], normVertCoord_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[2], xKernel_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[3], yKernel_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[4], zKernel_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[5], xNorm_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[6], yNorm_.data()))) ERR(retval);
  if ((retval = nc_put_var_double(varIds[0], varIds[7], zNorm_.data()))) ERR(retval);

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
    // No horizontal interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
          for (size_t jsv = 0; jsv < verStencilSize_[jz0]; ++jsv) {
            modelView(jnode0, jz0) += verWeights_[jz0][jsv]*redView(jnode0, verStencil_[jz0][jsv]);
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
      for (size_t jz = 0; jz < nz_; ++jz) {
        rSendVec[js*nz_+jz] = redView(jnode, jz);
      }
    }

    // Communication
    std::vector<double> mRecvVec(mRecvSize_*nz_);
    comm_.allToAllv(rSendVec.data(), rSendCounts3D.data(), rSendDispls3D.data(),
      mRecvVec.data(), mRecvCounts3D.data(), mRecvDispls3D.data());

    // Interpolation
    for (size_t jnode0 = 0; jnode0 < mSize_; ++jnode0) {
      for (size_t jsh = 0; jsh < horStencilSize_[jnode0]; ++jsh) {
        for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
          for (size_t jsv = 0; jsv < verStencilSize_[jz0]; ++jsv) {
            const size_t mIndex = horStencil_[jnode0][jsh]*nz_+verStencil_[jz0][jsv];
            modelView(jnode0, jz0) += horWeights_[jnode0][jsh]*verWeights_[jz0][jsv]
              *mRecvVec[mIndex];
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
        for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
          for (size_t jsv = 0; jsv < verStencilSize_[jz0]; ++jsv) {
            redView(jnode0, verStencil_[jz0][jsv]) += verWeights_[jz0][jsv]*modelView(jnode0, jz0);
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
      for (size_t jsh = 0; jsh < horStencilSize_[jnode0]; ++jsh) {
        for (size_t jz0 = 0; jz0 < nz0_; ++jz0) {
          for (size_t jsv = 0; jsv < verStencilSize_[jz0]; ++jsv) {
            const size_t mIndex = horStencil_[jnode0][jsh]*nz_+verStencil_[jz0][jsv];
            mRecvVec[mIndex] += horWeights_[jnode0][jsh]*verWeights_[jz0][jsv]
              *modelView(jnode0, jz0);
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
      for (size_t jz = 0; jz < nz_; ++jz) {
        redView(jnode, jz) += rSendVec[js*nz_+jz];
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
