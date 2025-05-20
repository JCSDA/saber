/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierTransform.h"

#include <algorithm>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/RandomField.h"

#include "saber/bifourier/BifourierUtilities.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

BifourierTransform::BifourierTransform(const oops::GeometryData & gdata,
                                       const std::string & gridUid,
                                       const oops::Variables & activeVars,
                                       const eckit::Configuration & params) :
    gdata_(gdata),
    comm_(gdata_.comm()),
    myrank_(comm_.rank()),
    params_(params),
    gridUid_(gridUid)
{
  oops::Log::trace() << classname() << "::BifourierTransform starting" << std::endl;

  // Check function space type
  ASSERT(gdata_.functionSpace().type() == "StructuredColumns");

  // Print active variables
  oops::Log::info() << "Info     : New Bifourier transform" << std::endl;
  oops::Log::info() << "Info     : - Active variable: " << activeVars << std::endl;

  // Get function space
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());

  // Get grid size
  nx_ = fs.grid().nx()[0];
  ny_ = fs.grid().ny();
  nodes_ = fs.size();
  oops::Log::test() << "- Regional grid size: " << nx_ << "x" << ny_ << std::endl;

  // Cell size
  dx_ = fs.grid().dx(0);
  dy_ = fs.grid().y(1) - fs.grid().y(0);
  oops::Log::test() << "- Cell sizes: " << dx_*1.0e-3 << " km x " << dy_*1.0e-3 << " km"
    << std::endl;

  // Mean latitude
  atlas::PointLonLat p1 = fs.grid().lonlat(0, 0);
  atlas::PointLonLat p2 = fs.grid().lonlat(nx_-1, ny_-1);
  meanLat_ = 0.5*(p1[1]+p2[1]);
  oops::Log::test() << "- Mean latitude: " << meanLat_ << " deg" << std::endl;

  // Number of levels for all variables
  nvz_ = 0;
  for (const auto & var : activeVars) {
    nvz_ += var.getLevels();
  }

  // Setup global spectral space parameters
  setupGlobalSpectralSpace();

  // Setup parallelization
  setupParallelization();

  // Setup local spectral space
  setupLocalSpectralSpace();

  // Setup FFT
  setupFFT();

  if (!params_.getBool("skip tests")) {
    // Get tests tolerance
    const double tolerance = params_.getDouble("spectral tolerance");

    // Generate random FieldSet
    atlas::FieldSet gpFset = util::createRandomFieldSet(comm_, gdata_.functionSpace(), activeVars);

    // Truncate grid-point field
    atlas::FieldSet spFset;
    gp2sp(gpFset, spFset, activeVars);
    sp2gp(spFset, gpFset, activeVars);

    // Grid-point to spectral
    gp2sp(gpFset, spFset, activeVars);

    // Check inverse
    atlas::FieldSet gpFsetTest = util::copyFieldSet(gpFset);
    sp2gp(spFset, gpFsetTest, activeVars);
    ASSERT(util::compareFieldSets(comm_, gpFset, gpFsetTest));
    oops::Log::test() << "- Direct-inverse test passed" << std::endl;

    // Check forward
    atlas::FieldSet spFsetTest;
    gp2sp(gpFsetTest, spFsetTest, activeVars);
    for (const auto & var : activeVars) {
      const size_t nz = var.getLevels();
      const auto spField = spFset[var.name()];
      const auto spFieldTest = spFsetTest[var.name()];
      const auto spView = make_view<double, 2>(spField);
      const auto spViewTest = make_view<double, 2>(spFieldTest);
      for (size_t js = 0; js < ns_; ++js) {
        for (size_t jz = 0; jz < nz; ++jz) {
          ASSERT(oops::is_close_relative(spView(js, jz), spViewTest(js, jz), tolerance));
        }
      }
    }
    oops::Log::test() << "- Inverse-direct test passed" << std::endl;

    // Check Parseval's identity
    double gpSqNorm = util::dotProductFieldSets(gpFset, gpFset, activeVars.variables(), comm_);
    double spSqNorm = util::dotProductFieldSets(spFset, spFset, activeVars.variables(), comm_);
    ASSERT(oops::is_close_relative(gpSqNorm, spSqNorm, tolerance));
    oops::Log::test() << "- Parseval identity test passed" << std::endl;

    // Adjoint test
    createRandomSpectralFieldSet(spFsetTest, activeVars);
    sp2gp(spFsetTest, gpFsetTest, activeVars);
    gpSqNorm = util::dotProductFieldSets(gpFset, gpFsetTest, activeVars.variables(), comm_);
    spSqNorm = util::dotProductFieldSets(spFset, spFsetTest, activeVars.variables(), comm_);
    ASSERT(oops::is_close_relative(gpSqNorm, spSqNorm, tolerance));
    oops::Log::test() << "- Adjoint test passed" << std::endl;

    // Derivatives / Laplacian consistency test
    for (const auto & var : activeVars) {
      // Get field
      const size_t nz = var.getLevels();
      auto spField = spFset[var.name()];

      // Double derivative in X direction
      atlas::Field spDxField;
      derivative(spField, spDxField, "x");
      atlas::Field spDx2Field;
      derivative(spDxField, spDx2Field, "x");

      // Double derivative in Y direction
      atlas::Field spDyField;
      derivative(spField, spDyField, "y");
      atlas::Field spDy2Field;
      derivative(spDyField, spDy2Field, "y");

      // Direct Laplacian
      atlas::Field spLapDirField = spField.clone();
      directLaplacian(spLapDirField);

      // Comparison
      const auto spDx2View = make_view<double, 2>(spDx2Field);
      const auto spDy2View = make_view<double, 2>(spDy2Field);
      const auto spLapDirView = make_view<double, 2>(spLapDirField);
      for (size_t js = 0; js < ns_; ++js) {
        for (size_t jz = 0; jz < nz; ++jz) {
          ASSERT(oops::is_close_relative(spLapDirView(js, jz),
            spDx2View(js, jz) + spDy2View(js, jz), tolerance));
        }
      }
    }
    oops::Log::test() << "- Derivatives / direct Laplacian consistency test passed" << std::endl;
  }

  oops::Log::trace() << classname() << "::BifourierTransform done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::gp2sp(const atlas::FieldSet & gpFset,
                               atlas::FieldSet & spFset,
                               const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::gp2sp starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create send and recv vectors
  std::vector<double> sendVec;
  std::vector<double> recvVec;

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Serialize from grid-point FieldSet
  recvVec.resize(gridRecvSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto gpField = gpFset[var.name()];
    ASSERT(gpField.shape(0) == static_cast<int>(nodes_));
    ASSERT(gpField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          recvVec[jgrv] = gpView(jnode, jz);
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  sendVec.resize(rowsSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    sendVec.data(), rowsSendCounts_.data(), rowsSendDispls_.data());

  // Reserialize
  for (size_t jrs = 0; jrs < rowsSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      size_t jf = rowsSendIndex_[jrs] + jvz*nx_;

      // Copy data
      rowsBufR_[jf] = sendVec[jrsv];
    }
  }

  // Compute direct transform
  fftw_execute(rowsPlan_r2c_);

  // Reserialize
  recvVec.resize(rowsRecvSize_*nvz_*2);
  size_t jrr = 0;
  for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
    for (size_t jk = 0; jk < nk_; ++jk) {
      for (size_t jvz = 0; jvz < nvz_; ++jvz) {
        // FFT vector index
        const size_t jf = jy*nvz_*nk_ + jvz*nk_ + jk;

        // Communication vector index
        const size_t jrrv = rowsRecvIndex_[jrr] + jvz*2;

        // Copy data
        recvVec[jrrv] = rowsBufC_[jf][0];
        recvVec[jrrv+1] = rowsBufC_[jf][1];
      }
      ++jrr;
    }
  }

  // Communication
  sendVec.resize(colsSendSize_*nvz_*2);
  comm_.allToAllv(recvVec.data(), rowsRecvCounts_.data(), rowsRecvDispls_.data(),
    sendVec.data(), colsSendCounts_.data(), colsSendDispls_.data());

  // Reserialize
  for (size_t jcs = 0; jcs < colsSendSize_; ++jcs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jcsv = jcs*nvz_*2 + jvz*2;

      // FFT vector index
      const size_t jf = colsSendIndex_[jcs] + jvz*2*ny_;

      // Copy data
      colsBufR_[jf] = sendVec[jcsv];
      colsBufR_[jf+ny_] = sendVec[jcsv+1];
    }
  }

  // Compute direct transform
  fftw_execute(colsPlan_r2c_);

  // Reserialize
  recvVec.resize(colsRecvSize_*nvz_);
  size_t jj = 0;
  for (size_t jk = 0; jk < nkPerTask_[myrank_]; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
       for (size_t jc = 0; jc < colsJq_[jk*nl_+jl].size(); ++jc) {
         // Get jq
         const size_t jq = colsJq_[jk*nl_+jl][jc];

         for (size_t jvz = 0; jvz < nvz_; ++jvz) {
          // FFT vector index
          const size_t jf = jk*nvz_*2*nl_ + jvz*2*nl_ + jl;

          // Communication vector index
          const size_t jcrv = colsRecvIndex_[jj] + jvz;

          // Copy data
          if (jq == 0) {
            recvVec[jcrv] = colsBufC_[jf][0];
          }
          if (jq == 1) {
            recvVec[jcrv] = colsBufC_[jf][1];
          }
          if (jq == 2) {
            recvVec[jcrv] = colsBufC_[jf+nl_][0];
          }
          if (jq == 3) {
            recvVec[jcrv] = colsBufC_[jf+nl_][1];
          }
        }

        // Update communication vector index
        ++jj;
      }
    }
  }

  // Communication
  sendVec.resize(eqchSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), colsRecvCounts_.data(), colsRecvDispls_.data(),
    sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data());

  // Prepare spectral FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (spFset.has(var.name())) {
      // Check sizes
      ASSERT(spFset[var.name()].shape(0) == static_cast<int>(ns_));
      ASSERT(spFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field spField = spFspace_->createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      spFset.add(spField);
    }
  }

  // Reserialize into spectral FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto spField = spFset[var.name()];

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        spView(js, jz) = sendVec[jesv];

        // Normalize FFT
        spView(js, jz) *= std::sqrt(normFFT_*spNorm(js));
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::gp2sp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::sp2gp(const atlas::FieldSet & spFset,
                               atlas::FieldSet & gpFset,
                               const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::sp2gp starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create send and recv vectors
  std::vector<double> sendVec;
  std::vector<double> recvVec;

  // Reserialize from spectral FieldSet
  sendVec.resize(eqchSendSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = spFset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        sendVec[jesv] = spView(js, jz);

        // Normalize FFT
        sendVec[jesv] *= std::sqrt(normFFT_/spNorm(js));
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  recvVec.resize(colsRecvSize_*nvz_);
  comm_.allToAllv(sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data(),
    recvVec.data(), colsRecvCounts_.data(), colsRecvDispls_.data());

  // Set FFT vector to zero
  for (size_t jj = 0; jj < nkPerTask_[myrank_]*nl_*nvz_*2; ++jj) {
    colsBufC_[jj][0] = 0.0;
    colsBufC_[jj][1] = 0.0;
  }

  // Reserialize
  size_t jj = 0;
  for (size_t jk = 0; jk < nkPerTask_[myrank_]; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
       for (size_t jc = 0; jc < colsJq_[jk*nl_+jl].size(); ++jc) {
         // Get jq
         const size_t jq = colsJq_[jk*nl_+jl][jc];

         for (size_t jvz = 0; jvz < nvz_; ++jvz) {
          // FFT vector index
          const size_t jf = jk*nvz_*2*nl_ + jvz*2*nl_ + jl;

          // Communication vector index
          const size_t jcrv = colsRecvIndex_[jj] + jvz;

          // Copy data
          if (jq == 0) {
            colsBufC_[jf][0] = recvVec[jcrv];
          }
          if (jq == 1) {
            colsBufC_[jf][1] = recvVec[jcrv];
          }
          if (jq == 2) {
            colsBufC_[jf+nl_][0] = recvVec[jcrv];
          }
          if (jq == 3) {
            colsBufC_[jf+nl_][1] = recvVec[jcrv];
          }
        }

        // Update communication vector index
        ++jj;
      }
    }
  }

  // Compute inverse transform
  fftw_execute(colsPlan_c2r_);

  // Reserialize
  sendVec.resize(colsSendSize_*nvz_*2);
  for (size_t jcs = 0; jcs < colsSendSize_; ++jcs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jcsv = jcs*nvz_*2 + jvz*2;

      // FFT vector index
      const size_t jf = colsSendIndex_[jcs] + jvz*2*ny_;

      // Copy data
      sendVec[jcsv] = colsBufR_[jf];
      sendVec[jcsv+1] = colsBufR_[jf+ny_];
    }
  }

  // Communication
  recvVec.resize(rowsRecvSize_*nvz_*2);
  comm_.allToAllv(sendVec.data(), colsSendCounts_.data(), colsSendDispls_.data(),
    recvVec.data(), rowsRecvCounts_.data(), rowsRecvDispls_.data());

  // Reserialize
  size_t jrr = 0;
  for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
    for (size_t jk = 0; jk < nk_; ++jk) {
      for (size_t jvz = 0; jvz < nvz_; ++jvz) {
        // FFT vector index
        const size_t jf = jy*nvz_*nk_ + jvz*nk_ + jk;

        // Communication vector index
        const size_t jrrv = rowsRecvIndex_[jrr] + jvz*2;

        // Copy data
        rowsBufC_[jf][0] = recvVec[jrrv];
        rowsBufC_[jf][1] = recvVec[jrrv+1];
      }
      ++jrr;
    }
  }

  // Compute inverse transform
  fftw_execute(rowsPlan_c2r_);

  // Reserialize
  sendVec.resize(rowsSendSize_*nvz_);
  for (size_t jrs = 0; jrs < rowsSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      const size_t jf = rowsSendIndex_[jrs] + jvz*nx_;

      // Copy data
      sendVec[jrsv] = rowsBufR_[jf];
    }
  }

  // Communication
  recvVec.resize(gridRecvSize_*nvz_);
  comm_.allToAllv(sendVec.data(), rowsSendCounts_.data(), rowsSendDispls_.data(),
    recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data());

  // Prepare grid-point FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (gpFset.has(var.name())) {
      // Check sizes
      ASSERT(gpFset[var.name()].shape(0) == static_cast<int>(nodes_));
      ASSERT(gpFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field gpField = gdata_.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      gpFset.add(gpField);
    }
  }

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Deserialize into grid-point FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto gpField = gpFset[var.name()];

    // Get field view
    auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          gpView(jnode, jz) = recvVec[jgrv];
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::sp2gp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::createRandomSpectralFieldSet(atlas::FieldSet & spFset,
                                                      const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::createRandomSpectralFieldSet starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Global vector
  std::vector<double> rand_vec_glb;

  if (myrank_ == 0) {
    // Generate global random vector
    rand_vec_glb.resize(nsGlb_*nvz_);
    util::NormalDistributionField dist(nsGlb_*nvz_, 0.0, 1.0);
    for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
      for (size_t jvz = 0; jvz < nvz_; ++jvz) {
        const size_t jj = jsGlb*nvz_ + jvz;
        const size_t jjOrdered = sMapping_[jsGlb]*nvz_ + jvz;
        rand_vec_glb[jj] = dist[jjOrdered];
      }
    }
  }

  // Scatter random vector
  std::vector<int> counts = sCounts_;
  std::vector<int> displs = sDispls_;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    counts[jt] *= nvz_;
    displs[jt] *= nvz_;
  }
  std::vector<double> rand_vec_loc(ns_*nvz_);
  comm_.scatterv(rand_vec_glb.cbegin(), rand_vec_glb.cend(), counts, displs,
    rand_vec_loc.begin(), rand_vec_loc.end(), 0);

  // Prepare spectral FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (spFset.has(var.name())) {
      // Check sizes
      ASSERT(spFset[var.name()].shape(0) == static_cast<int>(ns_));
      ASSERT(spFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field spField = spFspace_->createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      spFset.add(spField);
    }
  }

  // Reserialize into spectral FieldSet
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // CHeck field
    auto spField = spFset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t js = 0; js < ns_; ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Random vector index
        const size_t jr = js*nvz_ + jvz;

        // Copy data
        spView(js, jz) = rand_vec_loc[jr];
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::createRandomSpectralFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::derivative(const atlas::Field & spField,
                                    atlas::Field & spDerivField,
                                    const std::string & direction,
                                    const bool & adjoint) const {
  oops::Log::trace() << classname() << "::derivative starting" << std::endl;

  // Check field size
  ASSERT(spField.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = spField.shape(1);

  if (spDerivField.valid()) {
    // spDerivField is already allocated
    ASSERT(spDerivField.shape(0) == static_cast<int>(ns_));
    ASSERT(spDerivField.shape(1) == static_cast<int>(nz));
  } else {
    // Allocate spDerivField
    spDerivField = spField.clone();
  }

  // Get fields views
  const auto spView = make_view<double, 2>(spField);
  auto spDerivView = make_view<double, 2>(spDerivField);

  // Get derivative linear operator
  const auto firstIndexView = adjoint ?
    make_view<int, 1>(derivatives_[direction + "DerivCol"]) :
    make_view<int, 1>(derivatives_[direction + "DerivRow"]);
  const auto secondIndexView = adjoint ?
    make_view<int, 1>(derivatives_[direction + "DerivRow"]) :
    make_view<int, 1>(derivatives_[direction + "DerivCol"]);
  const auto factorView = make_view<double, 1>(derivatives_[direction + "DerivS"]);

  // Apply derivative linear operator
  spDerivView.assign(0.0);
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      spDerivView(firstIndexView(js), jz) = spView(secondIndexView(js), jz)*factorView(js);
    }
  }

  oops::Log::trace() << classname() << "::derivative done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::directLaplacian(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::directLaplacian starting" << std::endl;

  // Check field size
  ASSERT(field.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = field.shape(1);

  // Get field view
  auto view = make_view<double, 2>(field);

  // Apply direct Laplacian factor
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      view(js, jz) *= lapDirVec_[js];
    }
  }

  oops::Log::trace() << classname() << "::directLaplacian done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::inverseLaplacian(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::inverseLaplacian starting" << std::endl;

  // Check field size
  ASSERT(field.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = field.shape(1);

  // Get field view
  auto view = make_view<double, 2>(field);

  // Apply inverse Laplacian factor
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      view(js, jz) *= lapInvVec_[js];
    }
  }
  oops::Log::trace() << classname() << "::inverseLaplacian done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::cv2fset(const atlas::Field & cv,
                                 atlas::FieldSet & fset,
                                 const oops::Variables & activeVars,
                                 const size_t & offset) const {
  oops::Log::trace() << classname() << "::cv2fset starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Clear FieldSet
  fset.clear();

  // Get control vector view
  const auto cvView = make_view<double, 1>(cv);

  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Create field
    atlas::Field spField = spFspace_->createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    fset.add(spField);

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t jz = 0; jz < nz; ++jz) {
      // Total level index
      const size_t jvz = zOffset + jz;

      for (size_t js = 0; js < ns_; ++js) {
        // Control vector index
        const size_t jcv = offset + jvz*ns_ + js;

        // Copy data
        spView(js, jz) = cvView(jcv);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::cv2fset done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::fset2cv(const atlas::FieldSet & fset,
                                 atlas::Field & cv,
                                 const oops::Variables & activeVars,
                                 const size_t & offset) const {
  oops::Log::trace() << classname() << "::fset2cv starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Get CV view
  auto cvView = make_view<double, 1>(cv);

  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = fset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jz = 0; jz < nz; ++jz) {
      // Total level index
      const size_t jvz = zOffset + jz;

      for (size_t js = 0; js < ns_; ++js) {
        // Control vector index
        const size_t jcv = offset + jvz*ns_ + js;

        // Copy data
        cvView(jcv) = spView(js, jz);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::fset2cv done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::copyFieldSet(const atlas::FieldSet & spInputFset,
                                      atlas::FieldSet & spOutputFset,
                                      const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::copyFieldSet starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Remove active variables from output FieldSet
  util::removeFieldsFromFieldSet(spOutputFset, activeVars.variables());

  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check input field
    const auto spInputField = spInputFset[var.name()];
    ASSERT(spInputField.shape(0) == static_cast<int>(ns_));
    ASSERT(spInputField.shape(1) == static_cast<int>(nz));

    // Create output field
    atlas::Field spOutputField = spFspace_->createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    spOutputFset.add(spOutputField);

    // Get fields views
    const auto spInputView = make_view<double, 2>(spInputField);
    auto spOutputView = make_view<double, 2>(spOutputField);

    for (size_t jz = 0; jz < nz; ++jz) {
      for (size_t js = 0; js < ns_; ++js) {
        // Copy data
        spOutputView(js, jz) = spInputView(js, jz);
      }
    }
  }

  oops::Log::trace() << classname() << "::copyFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

double BifourierTransform::kstar(const size_t & jk,
                                 const size_t & jl,
                                 const size_t & nk,
                                 const size_t & nl,
                                 const size_t & nwGlb) const {
  const double w = static_cast<double>(nwGlb-1)*std::sqrt(
    static_cast<double>(jk*jk)/static_cast<double>((nk-2)*(nk-2))
    + static_cast<double>(jl*jl)/static_cast<double>((nl-2)*(nl-2)));
  return w;
}

// -----------------------------------------------------------------------------

void BifourierTransform::print(std::ostream & os) const {
  os << classname();
}


// -----------------------------------------------------------------------------

void BifourierTransform::setupGlobalSpectralSpace() {
  oops::Log::trace() << classname() << "::setupGlobalSpectralSpace starting" << std::endl;

  // Maximum wave numbers
  nk_ = nx_/2+1;
  nl_ = ny_/2+1;
  oops::Log::test() << "- Spectral sizes: " << nk_ << "x" << nl_ << std::endl;

  // Differential operators factors
  exwn_ = 2.0*M_PI/(static_cast<double>(nx_)*dx_);
  eywn_ = 2.0*M_PI/(static_cast<double>(ny_)*dy_);

  // Normalization factor
  normFFT_ = 1.0/static_cast<double>(nx_*ny_);

  // Define truncation
  if (params_.getString("truncation type") == "arome") {
    // Same as the AROME model

    // Define tolerance to define jwGlb
    jwGlbTol_ = 0.49;

    // Define truncation parameters
    oops::Log::test() << "- Truncation parameters MxN: " << nk_-2 << "x" <<  nl_-2 << std::endl;

    // Define ellips
    ellips_.resize(nk_-1);
    ellips_[0] = nl_-2;
    for (size_t jk = 1; jk < nk_-2; ++jk) {
      ellips_[jk] = static_cast<size_t>(static_cast<float>(nl_-2)/static_cast<float>(nk_-2)
        *std::sqrt(static_cast<float>((nk_-2)*(nk_-2)-jk*jk))+1.0e-10);
    }
    ellips_[nk_-2] = 0;

    // Maximum total wave number
    nwGlb_ = std::max(nk_-1, nl_-1);
    oops::Log::test() << "- Maximum total wave number: " << nwGlb_-1 << std::endl;
  } else {
    // Unknown truncation
    throw eckit::Exception("unknown truncation type", Here());
  }

  // Mapping
  spNormKL_.resize(nk_*nl_);
  std::fill(spNormKL_.begin(), spNormKL_.end(), 0.0);
  size_t jk;
  size_t jl;

  // k = 0
  jk = 0;

  // l = 0
  jl = 0;
  addSpectralCoefficient(jk, jl, ReRe, 0, 0);

  // 0 < l < nl_-1
  for (size_t jl = 1; jl < nl_-1; ++jl) {
    addSpectralCoefficient(jk, jl, ReRe, 0, 1);
    addSpectralCoefficient(jk, jl, ReIm, 0, -1);
  }

  // l = nl_-1
  jl = nl_-1;
  if (ny_ % 2 == 1) {
    addSpectralCoefficient(jk, jl, ReRe, 0, 1);
    addSpectralCoefficient(jk, jl, ReIm, 0, -1);
  } else {
    addSpectralCoefficient(jk, jl, ReRe, 0, 0);
  }

  // 0 < k < nk_-1
  for (size_t jk = 1; jk < nk_-1; ++jk) {
    // l = 0
    jl = 0;
    addSpectralCoefficient(jk, jl, ReRe, 1, 0);
    addSpectralCoefficient(jk, jl, ImRe, -1, 0);

    // 0 < l < nl_-1
    for (size_t jl = 1; jl < nl_-1; ++jl) {
      addSpectralCoefficient(jk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jk, jl, ImIm, -2, -1);
    }

    // l = nl_-1
    jl = nl_-1;
    if (ny_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jk, jl, ImIm, -2, -1);
    } else {
      addSpectralCoefficient(jk, jl, ReRe, 1, 0);
      addSpectralCoefficient(jk, jl, ImRe, -1, 0);
    }
  }

  // k = nk_-1
  jk = nk_-1;

  // l = 0
  jl = 0;
  if (nx_ % 2 == 1) {
    addSpectralCoefficient(jk, jl, ReRe, 1, 0);
    addSpectralCoefficient(jk, jl, ImRe, -1, 0);
  } else {
    addSpectralCoefficient(jk, jl, ReRe, 0, 0);
  }

  // 0 < l < nl_-1
  for (size_t jl = 1; jl < nl_-1; ++jl) {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jk, jl, ImIm, -2, -1);
    } else {
      addSpectralCoefficient(jk, jl, ReRe, 0, 1);
      addSpectralCoefficient(jk, jl, ReIm, 0, -1);
    }
  }

  // l = nl_-1
  jl = nl_-1;
  if (ny_ % 2 == 1) {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jk, jl, ImIm, -2, -1);

    } else {
      addSpectralCoefficient(jk, jl, ReRe, 0, 1);
      addSpectralCoefficient(jk, jl, ReIm, 0, -1);
    }
  } else {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 1, 0);
      addSpectralCoefficient(jk, jl, ImRe, -1, 0);
    } else {
      addSpectralCoefficient(jk, jl, ReRe, 0, 0);
    }
  }

  // Define global spectral size
  nsGlb_ = spVec_.size();
  oops::Log::test() << "- Spectral array global size: " << nsGlb_ << std::endl;

  oops::Log::trace() << classname() << "::setupGlobalSpectralSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::setupParallelization() {
  oops::Log::trace() << classname() << "::setupParallelization starting" << std::endl;

  // Split in y direction
  nyPerTask_.resize(comm_.size());
  std::fill(nyPerTask_.begin(), nyPerTask_.end(), 0);
  size_t index = 0;
  for (size_t jy = 0; jy < ny_; ++jy) {
    ++nyPerTask_[index];
    ++index;
    if (index == comm_.size()) index = 0;
  }
  std::vector<size_t> nyStart(comm_.size());
  std::vector<size_t> nyEnd(comm_.size());
  nyStart[0] = 0;
  nyEnd[0] = nyPerTask_[0]-1;
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    nyStart[jt+1] = nyStart[jt]+nyPerTask_[jt];
    nyEnd[jt+1] = nyStart[jt+1]+nyPerTask_[jt+1]-1;
  }

  // Split in k direction
  nkPerTask_.resize(comm_.size());
  std::fill(nkPerTask_.begin(), nkPerTask_.end(), 0);
  index = 0;
  for (size_t jk = 0; jk < nk_; ++jk) {
    ++nkPerTask_[index];
    ++index;
    if (index == comm_.size()) index = 0;
  }
  std::vector<size_t> nkStart(comm_.size());
  std::vector<size_t> nkEnd(comm_.size());
  nkStart[0] = 0;
  nkEnd[0] = nkPerTask_[0]-1;
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    nkStart[jt+1] = nkStart[jt]+nkPerTask_[jt];
    nkEnd[jt+1] = nkStart[jt+1]+nkPerTask_[jt+1]-1;
  }

  // Split truncation wavenumbers in equal chunks
  std::vector<int> jkVec;
  std::vector<int> jlVec;
  std::vector<int> jwGlbVec;
  std::vector<size_t> nklPerTaskTarget(comm_.size());
  std::fill(nklPerTaskTarget.begin(), nklPerTaskTarget.end(), 0);
  index = 0;
  for (size_t jk = 0; jk < ellips_.size(); ++jk) {
    for (size_t jl = 0; jl <= ellips_[jk]; ++jl) {
      jkVec.push_back(jk);
      jlVec.push_back(jl);
      jwGlbVec.push_back(kstar(jk, jl, nk_, nl_, nwGlb_));
      ++nklPerTaskTarget[index];
      ++index;
      if (index == comm_.size()) index = 0;
    }
  }

  // Sort truncation wavenumbers in ascending jwGlb
  const size_t nkl = jwGlbVec.size();
  std::vector<size_t> klOrder(nkl);
  std::iota(klOrder.begin(), klOrder.end(), 0);
  std::stable_sort(klOrder.begin(), klOrder.end(),
    [&](int i, int j){return jwGlbVec[i] < jwGlbVec[j];});

  // Split total wave number among tasks
  std::vector<size_t> nklPerTask(comm_.size());
  std::vector<std::vector<size_t>> ellipsTask(ellips_.size());
  for (size_t jk = 0; jk < ellips_.size(); ++jk) {
    ellipsTask[jk].resize(ellips_[jk]+1);
  }
  size_t jt = 0;
  for (size_t jkl = 0; jkl < nkl; ++jkl) {
    // Conditions to switch to the next task
    if ((nklPerTask[jt] > nklPerTaskTarget[jt]) && (jt < comm_.size()-1)) {
      nklPerTaskTarget[jt+1] += nklPerTaskTarget[jt] - nklPerTask[jt];
      ++jt;
    }

    // Get jk, jl
    const size_t jk = jkVec[klOrder[jkl]];
    const size_t jl = jlVec[klOrder[jkl]];

    // Save task
    ellipsTask[jk][jl] = jt;

    // Update local size
    ++nklPerTask[jt];
  }

  // Add task in global spectral vector
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get jk, jl
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;

    // Add task
    spVec_[jsGlb].jt = ellipsTask[jk][jl];
  }

  // Spectral size per task and local to global mapping
  nsPerTask_.resize(comm_.size());
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get task
    const size_t jt = spVec_[jsGlb].jt;

    // Update number of spectral coefficient for this task
    ++nsPerTask_[jt];

    if (jt == myrank_) {
      // Add local spectral coefficient
      sToSGlb_.push_back(jsGlb);
    }
  }

  // Save local size
  ns_ = sToSGlb_.size();

  // Communication vectors and mapping
  sCounts_.resize(comm_.size());
  sDispls_.resize(comm_.size());
  if (myrank_ == 0) {
    sMapping_.resize(nsGlb_);
  }
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sCounts_[jt] = nsPerTask_[jt];
    sDispls_[jt] = static_cast<int>(jt ? sDispls_[jt-1] + sCounts_[jt-1] : 0);
  }
  comm_.gatherv(sToSGlb_.cbegin(), sToSGlb_.cend(), sMapping_.begin(), sMapping_.end(),
    sCounts_, sDispls_, 0);


  // Compute spectral imbalance
  const double sImb = static_cast<double>(*std::max_element(nsPerTask_.begin(), nsPerTask_.end()))
    / static_cast<double>(*std::min_element(nsPerTask_.begin(), nsPerTask_.end()));
  oops::Log::info() << "Info     : - Spectral imbalance (max/min): " << sImb << std::endl;

  // Rows <=> grid

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Index fields views
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  const auto indexIView = make_view<int, 1>(fs.index_i());
  const auto indexJView = make_view<int, 1>(fs.index_j());

  // Number of values on each destination task
  gridRecvSize_ = 0;
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    if (ghostView(jnode) == 0) {
      ++gridRecvSize_;
    }
  }

  // Define destination task
  std::vector<int> rowsSendTask(gridRecvSize_);
  std::vector<int> rowsSendOffset(gridRecvSize_);
  std::vector<int> rowsSendOffsetPerTask(comm_.size(), 0);
  size_t jgr = 0;
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    if (ghostView(jnode) == 0) {
      for (size_t jt = 0; jt < comm_.size(); ++jt) {
        if (static_cast<size_t>(indexJView(jnode))-1 >= nyStart[jt] &&
          static_cast<size_t>(indexJView(jnode))-1 <= nyEnd[jt]) {
          rowsSendTask[jgr] = jt;
          rowsSendOffset[jgr] = rowsSendOffsetPerTask[jt];
          ++rowsSendOffsetPerTask[jt];
        }
      }
      ++jgr;
    }
  }

  // RecvCounts
  gridRecvCounts_.resize(comm_.size());
  std::fill(gridRecvCounts_.begin(), gridRecvCounts_.end(), 0);
  for (size_t jgr = 0; jgr < gridRecvSize_; ++jgr) {
    ++gridRecvCounts_[rowsSendTask[jgr]];
  }

  // RecvDispls
  gridRecvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    gridRecvDispls_[jt] = static_cast<int>(jt ? gridRecvDispls_[jt-1] + gridRecvCounts_[jt-1] : 0);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> rRecvCountsBuffer(comm_.size());
  comm_.allGatherv(gridRecvCounts_.begin(), gridRecvCounts_.end(), rRecvCountsBuffer);
  std::vector<int> rRecvCountsGlb_ = std::move(rRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rowsSendCounts_.push_back(rRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  rowsSendSize_ = 0;
  for (const auto & n : rowsSendCounts_) rowsSendSize_ += n;

  // SendDispls
  rowsSendDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rowsSendDispls_[jt] = static_cast<int>(jt ? rowsSendDispls_[jt-1] + rowsSendCounts_[jt-1] : 0);
  }

  // Communicate indices
  std::vector<int> gridRecvIndex_x(gridRecvSize_);
  std::vector<int> gridRecvIndex_y(gridRecvSize_);
  jgr = 0;
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    if (ghostView(jnode) == 0) {
      gridRecvIndex_x[jgr] = indexIView(jnode)-1;
      gridRecvIndex_y[jgr] = indexJView(jnode)-1;
      ++jgr;
    }
  }
  std::vector<int> rowsSendIndex_x(rowsSendSize_);
  std::vector<int> rowsSendIndex_y(rowsSendSize_);
  rowsSendIndex_.resize(rowsSendSize_);
  comm_.allToAllv(gridRecvIndex_x.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    rowsSendIndex_x.data(), rowsSendCounts_.data(), rowsSendDispls_.data());
  comm_.allToAllv(gridRecvIndex_y.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    rowsSendIndex_y.data(), rowsSendCounts_.data(), rowsSendDispls_.data());
  for (size_t jrs = 0; jrs < rowsSendSize_; ++jrs) {
    rowsSendIndex_[jrs] = (rowsSendIndex_y[jrs]-nyStart[myrank_])*nvz_*nx_ + rowsSendIndex_x[jrs];
  }

  // Effective index
  gridRecvIndex_.resize(gridRecvSize_);
  for (size_t jgr = 0; jgr < gridRecvSize_; ++jgr) {
    gridRecvIndex_[jgr] = (gridRecvDispls_[rowsSendTask[jgr]] + rowsSendOffset[jgr])*nvz_;
  }

  // Columns <=> rows

  // Number of values on each destination task
  rowsRecvSize_ = nyPerTask_[myrank_]*nk_;

  // Define destination task
  std::vector<int> colsSendTask(rowsRecvSize_);
  std::vector<int> colsSendOffset(rowsRecvSize_);
  std::vector<int> colsSendOffsetPerTask(comm_.size(), 0);
  for (size_t jk = 0; jk < nk_; ++jk) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      for (size_t jt = 0; jt < comm_.size(); ++jt) {
        if (jk >= nkStart[jt] && jk <= nkEnd[jt]) {
          colsSendTask[jy*nk_+jk] = jt;
          colsSendOffset[jy*nk_+jk] = colsSendOffsetPerTask[jt];
          ++colsSendOffsetPerTask[jt];
        }
      }
    }
  }

  // RecvCounts
  rowsRecvCounts_.resize(comm_.size());
  std::fill(rowsRecvCounts_.begin(), rowsRecvCounts_.end(), 0);
  for (size_t jrr = 0; jrr < rowsRecvSize_; ++jrr) {
    ++rowsRecvCounts_[colsSendTask[jrr]];
  }

  // RecvDispls
  rowsRecvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    rowsRecvDispls_[jt] = static_cast<int>(jt ? rowsRecvDispls_[jt-1] + rowsRecvCounts_[jt-1] : 0);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> rowsRecvCountsBuffer(comm_.size());
  comm_.allGatherv(rowsRecvCounts_.begin(), rowsRecvCounts_.end(), rowsRecvCountsBuffer);
  std::vector<int> rowsRecvCountsGlb_ = std::move(rowsRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    colsSendCounts_.push_back(rowsRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  colsSendSize_ = 0;
  for (const auto & n : colsSendCounts_) colsSendSize_ += n;

  // SendDispls
  colsSendDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    colsSendDispls_[jt] = static_cast<int>(jt ? colsSendDispls_[jt-1] + colsSendCounts_[jt-1] : 0);
  }

  // Communicate indices
  std::vector<int> rowsRecvIndex_k(rowsRecvSize_);
  std::vector<int> rowsRecvIndex_y(rowsRecvSize_);
  for (size_t jk = 0; jk < nk_; ++jk) {
    for (size_t jy = 0; jy < nyPerTask_[myrank_]; ++jy) {
      rowsRecvIndex_k[jk*nyPerTask_[myrank_]+jy] = jk;
      rowsRecvIndex_y[jk*nyPerTask_[myrank_]+jy] = jy+nyStart[myrank_];
    }
  }
  std::vector<int> colsSendIndex_k(colsSendSize_);
  std::vector<int> colsSendIndex_y(colsSendSize_);
  colsSendIndex_.resize(colsSendSize_);
  comm_.allToAllv(rowsRecvIndex_k.data(), rowsRecvCounts_.data(), rowsRecvDispls_.data(),
    colsSendIndex_k.data(), colsSendCounts_.data(), colsSendDispls_.data());
  comm_.allToAllv(rowsRecvIndex_y.data(), rowsRecvCounts_.data(), rowsRecvDispls_.data(),
    colsSendIndex_y.data(), colsSendCounts_.data(), colsSendDispls_.data());
  for (size_t jcs = 0; jcs < colsSendSize_; ++jcs) {
    colsSendIndex_[jcs] = (colsSendIndex_k[jcs]-nkStart[myrank_])*nvz_*2*ny_ + colsSendIndex_y[jcs];
  }

  // Effective index
  rowsRecvIndex_.resize(rowsRecvSize_);
  for (size_t jrr = 0; jrr < rowsRecvSize_; ++jrr) {
    rowsRecvIndex_[jrr] = (rowsRecvDispls_[colsSendTask[jrr]] + colsSendOffset[jrr])*nvz_*2;
  }

  // Rows <=> equal chunks

  // Number of values on each destination task
  colsRecvSize_ = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    if (jk >= nkStart[myrank_] && jk <= nkEnd[myrank_]) {
      ++colsRecvSize_;
    }
  }

  // Define destination task
  std::vector<int> eqchSendTask(colsRecvSize_);
  std::vector<int> eqchSendOffset(colsRecvSize_);
  std::vector<int> eqchSendOffsetPerTask(comm_.size(), 0);
  size_t jv = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    if (jk >= nkStart[myrank_] && jk <= nkEnd[myrank_]) {
      const size_t jt = spVec_[jsGlb].jt;
      eqchSendTask[jv] = jt;
      eqchSendOffset[jv] = eqchSendOffsetPerTask[jt];
      ++eqchSendOffsetPerTask[jt];
      ++jv;
    }
  }

  // RecvCounts
  colsRecvCounts_.resize(comm_.size());
  std::fill(colsRecvCounts_.begin(), colsRecvCounts_.end(), 0);
  for (size_t jcr = 0; jcr < colsRecvSize_; ++jcr) {
    ++colsRecvCounts_[eqchSendTask[jcr]];
  }

  // RecvDispls
  colsRecvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    colsRecvDispls_[jt] = static_cast<int>(jt ? colsRecvDispls_[jt-1] + colsRecvCounts_[jt-1] : 0);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> colsRecvCountsBuffer(comm_.size());
  comm_.allGatherv(colsRecvCounts_.begin(), colsRecvCounts_.end(), colsRecvCountsBuffer);
  std::vector<int> colsRecvCountsGlb_ = std::move(colsRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    eqchSendCounts_.push_back(colsRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  eqchSendSize_ = 0;
  for (const auto & n : eqchSendCounts_) eqchSendSize_ += n;

  // SendDispls
  eqchSendDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    eqchSendDispls_[jt] = static_cast<int>(jt ? eqchSendDispls_[jt-1] + eqchSendCounts_[jt-1] : 0);
  }

  // Get destination task inverse order
  std::vector<size_t> colsRecvOrder(colsRecvSize_);
  std::iota(colsRecvOrder.begin(), colsRecvOrder.end(), 0);
  std::stable_sort(colsRecvOrder.begin(), colsRecvOrder.end(),
    [&](int i, int j){return eqchSendTask[i] < eqchSendTask[j];});
  std::vector<size_t> colsRecvOrderInverse(colsRecvSize_);
  for (size_t jcr = 0; jcr < colsRecvSize_; ++jcr) {
    colsRecvOrderInverse[colsRecvOrder[jcr]] = jcr;
  }

  // Communicate indices
  std::vector<int> colsRecvIndex(colsRecvSize_);
  size_t jcr = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    if (jk >= nkStart[myrank_] && jk <= nkEnd[myrank_]) {
      colsRecvIndex[colsRecvOrderInverse[jcr]] = jsGlb;
      ++jcr;
    }
  }
  eqchSendIndex_.resize(eqchSendSize_);
  comm_.allToAllv(colsRecvIndex.data(), colsRecvCounts_.data(), colsRecvDispls_.data(),
    eqchSendIndex_.data(), eqchSendCounts_.data(), eqchSendDispls_.data());
  for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
    const auto it = std::find(sToSGlb_.begin(), sToSGlb_.end(), eqchSendIndex_[jes]);
    ASSERT(it != sToSGlb_.end());
    eqchSendIndex_[jes] = it-sToSGlb_.begin();
  }

  // Effective index
  colsRecvIndex_.resize(colsRecvSize_);
  for (size_t jcr = 0; jcr < colsRecvSize_; ++jcr) {
    colsRecvIndex_[jcr] = (colsRecvDispls_[eqchSendTask[jcr]] + eqchSendOffset[jcr])*nvz_;
  }

  // Scale counts and displs for all levels
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    gridRecvCounts_[jt] *= nvz_;
    gridRecvDispls_[jt] *= nvz_;
    rowsSendCounts_[jt] *= nvz_;
    rowsSendDispls_[jt] *= nvz_;
    rowsRecvCounts_[jt] *= nvz_*2;
    rowsRecvDispls_[jt] *= nvz_*2;
    colsSendCounts_[jt] *= nvz_*2;
    colsSendDispls_[jt] *= nvz_*2;
    colsRecvCounts_[jt] *= nvz_;
    colsRecvDispls_[jt] *= nvz_;
    eqchSendCounts_[jt] *= nvz_;
    eqchSendDispls_[jt] *= nvz_;
  }

  // Number of component for each local jk,jl couple
  colsJq_.resize(nkPerTask_[myrank_]*nl_);
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    if (jk >= nkStart[myrank_] && jk <= nkEnd[myrank_]) {
      const size_t jkPerTask = jk-nkStart[myrank_];
      const size_t jl = spVec_[jsGlb].jl;
      const size_t jq = spVec_[jsGlb].jq;
      colsJq_[jkPerTask*nl_+jl].push_back(jq);
    }
  }

  // Get min and max jwGlb
  nwStartPerTask_.resize(comm_.size(), nwGlb_-1);
  nwEndPerTask_.resize(comm_.size(), 0);
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get jwGlb and jt
    const size_t jwGlb = spVec_[jsGlb].jwGlb;
    const size_t jt = spVec_[jsGlb].jt;

    // Update min and max jwGlb
    nwStartPerTask_[jt] = std::min(nwStartPerTask_[jt], jwGlb);
    nwEndPerTask_[jt] = std::max(nwEndPerTask_[jt], jwGlb);
  }
  nwPerTask_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    nwPerTask_[jt] = nwEndPerTask_[jt] - nwStartPerTask_[jt] + 1;
  }

  // Local number of total wavenumbers
  nw_ = nwPerTask_[myrank_];

  // Define root nw
  nwRootPerTask_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    nwRootPerTask_[jt] = nwStartPerTask_[jt+1] - nwStartPerTask_[jt];
  }
  nwRootPerTask_[comm_.size()-1] = nwGlb_ - nwStartPerTask_[jt];

  // Local number of root total wavenumbers
  nwRoot_ = nwRootPerTask_[myrank_];

  // Communication vectors
  wCounts_.resize(comm_.size());
  wDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    wCounts_[jt] = nwRootPerTask_[jt];
    wDispls_[jt] = static_cast<int>(jt ? wDispls_[jt-1] + wCounts_[jt-1] : 0);
  }

  // Compute total wavenumber imbalance
  const double wImb = static_cast<double>(*std::max_element(nwPerTask_.begin(), nwPerTask_.end()))
    / static_cast<double>(*std::min_element(nwPerTask_.begin(), nwPerTask_.end()));
  oops::Log::info() << "Info     : - Total wavenumber imbalance (max/min): " << wImb << std::endl;

  oops::Log::trace() << classname() << "::setupParallelization done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::setupLocalSpectralSpace() {
  oops::Log::trace() << classname() << "::setupLocalSpectralSpace starting" << std::endl;

  // Create dummy PointCloud function space
  atlas::FieldSet flds;
  atlas::Field lonlatField("lonlat", make_datatype<double>(), make_shape(ns_, 2));
  atlas::Field ghostField("ghost", make_datatype<int>(), make_shape(ns_));
  atlas::Field remoteIndexField("remote_index", make_datatype<int>(), make_shape(ns_));
  atlas::Field partitionField("partition", make_datatype<int>(), make_shape(ns_));
  atlas::Field globalIndexField("global_index", make_datatype<int64_t>(), make_shape(ns_));
  flds.add(lonlatField);
  flds.add(ghostField);
  flds.add(remoteIndexField);
  flds.add(partitionField);
  flds.add(globalIndexField);
  auto lonlatView = make_view<double, 2>(lonlatField);
  auto ghostView = make_view<int, 1>(ghostField);
  auto remoteIndexView = make_view<int, 1>(remoteIndexField);
  auto partitionView = make_view<int, 1>(partitionField);
  auto globalIndexView = make_view<int64_t, 1>(globalIndexField);

  for (size_t js = 0; js < ns_; ++js) {
    // Get global index
    const size_t jsGlb = sToSGlb_[js];

    // Get jk/jl/jq
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;
    const size_t jq = spVec_[jsGlb].jq;

    // Define dummy lon/lat from jk/jl/jq
    lonlatView(js, 0) = (static_cast<double>(jk)+0.25*static_cast<double>(jq))
      /static_cast<double>(nk_)*360.0;
    lonlatView(js, 1) = -90.0*(static_cast<double>(jl)+0.25*static_cast<double>(jq))
      /static_cast<double>(nl_)*180.0;

    // Define ghost field
    ghostView(js) = 0;

    // Define remote index field
    remoteIndexView(js) = js;

    // Define partition field
    partitionView(js) = myrank_;

    // Define global index field
    globalIndexView(js) = static_cast<int64_t>(jsGlb);
  }
  spFspace_.reset(new atlas::functionspace::PointCloud(flds));

  // Generate spectral UID
  specUid_ = generateSpectralUid(*spFspace_, comm_);

  // Print UIDs
  oops::Log::info() << "Info     : - UIDs: " << gridUid_ << " / " << specUid_ << std::endl;

  // Allocate vectors
  jkVec_.resize(ns_);
  jlVec_.resize(ns_);
  jqVec_.resize(ns_);
  jwGlbVec_.resize(ns_);
  spNormVec_.resize(ns_);
  lapDirVec_.resize(ns_);
  lapInvVec_.resize(ns_);

  // Allocate fields
  atlas::Field xDerivRowField("xDerivRow", make_datatype<int>(), make_shape(ns_));
  atlas::Field xDerivColField("xDerivCol", make_datatype<int>(), make_shape(ns_));
  atlas::Field xDerivSField("xDerivS", make_datatype<double>(), make_shape(ns_));
  atlas::Field yDerivRowField("yDerivRow", make_datatype<int>(), make_shape(ns_));
  atlas::Field yDerivColField("yDerivCol", make_datatype<int>(), make_shape(ns_));
  atlas::Field yDerivSField("yDerivS", make_datatype<double>(), make_shape(ns_));

  // Add fields to derivatives_
  derivatives_.add(xDerivRowField);
  derivatives_.add(xDerivColField);
  derivatives_.add(xDerivSField);
  derivatives_.add(yDerivRowField);
  derivatives_.add(yDerivColField);
  derivatives_.add(yDerivSField);

  // Get fields views
  auto xDerivRowView = make_view<int, 1>(xDerivRowField);
  auto xDerivColView = make_view<int, 1>(xDerivColField);
  auto xDerivSView = make_view<double, 1>(xDerivSField);
  auto yDerivRowView = make_view<int, 1>(yDerivRowField);
  auto yDerivColView = make_view<int, 1>(yDerivColField);
  auto yDerivSView = make_view<double, 1>(yDerivSField);

  // Fill data
  for (size_t js = 0; js < ns_; ++js) {
    // Get global index
    const size_t jsGlb = sToSGlb_[js];

    // Get spVec_ values
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;
    const size_t jq = spVec_[jsGlb].jq;
    const size_t jwGlb = spVec_[jsGlb].jwGlb;
    const int jsXDerivativeOffset = spVec_[jsGlb].jsXDerivativeOffset;
    const int jsYDerivativeOffset = spVec_[jsGlb].jsYDerivativeOffset;

    // Copy jk, jl, jq and jwGlb
    jkVec_[js] = jk;
    jlVec_[js] = jl;
    jqVec_[js] = jq;
    jwGlbVec_[js] = jwGlb;

    // Spectral norm
    spNormVec_[js] = static_cast<double>(spNormKL_[jk*nl_+jl]);

    // X-derivative linear operator
    xDerivColView(js) = js;
    if (jsXDerivativeOffset == 0) {
      // Set to zero
      xDerivRowView(js) = js;
      xDerivSView(js) = 0.0;
    } else {
      // Involved in the derivative
      xDerivRowView(js) = js+jsXDerivativeOffset;
      if (jsXDerivativeOffset > 0) {
        xDerivSView(js) = static_cast<double>(jk)*exwn_;
      } else  {
        xDerivSView(js) = -static_cast<double>(jk)*exwn_;
      }
    }

    // Y-derivative linear operator
    yDerivColView(js) = js;
    if (jsYDerivativeOffset == 0) {
      // Set to zero
      yDerivRowView(js) = js;
      yDerivSView(js) = 0.0;
    } else {
      // Involved in the derivative
      yDerivRowView(js) = js+jsYDerivativeOffset;
      if (jsYDerivativeOffset > 0) {
        yDerivSView(js) = static_cast<double>(jl)*eywn_;
      } else  {
        yDerivSView(js) = -static_cast<double>(jl)*eywn_;
      }
    }

    // Direct Laplacian
    lapDirVec_[js] = -(static_cast<double>(jk*jk)*exwn_*exwn_
      + static_cast<double>(jl*jl)*eywn_*eywn_);

    // Inverse Laplacian
    if (jk == 0 && jl == 0) {
      lapInvVec_[js] = 0.0;
    } else {
      lapInvVec_[js] = 1.0/lapDirVec_[js];
    }
  }

  // Compute control vector size
  ctlVecSize_ = ns_*nvz_;

  oops::Log::trace() << classname() << "::setupLocalSpectralSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::setupFFT() {
  oops::Log::trace() << classname() << "::setupFFT starting" << std::endl;

  // Rows setup
  int xRank = 1;
  int xN[] = {static_cast<int>(nx_)};
  int xHowmany = nyPerTask_[myrank_]*nvz_;
  int *xInembed = NULL;
  const int xIstride = 1;
  const int xIdist = static_cast<int>(nx_);
  int *xOnembed = NULL;
  const int xOstride = 1;
  const int xOdist = static_cast<int>(nk_);
  rowsBufR_ = fftw_alloc_real(xHowmany*nx_);
  rowsBufC_ = fftw_alloc_complex(xHowmany*nk_);
  rowsPlan_r2c_ = fftw_plan_many_dft_r2c(xRank, xN, xHowmany, rowsBufR_, xInembed, xIstride, xIdist,
    rowsBufC_, xOnembed, xOstride, xOdist, FFTW_ESTIMATE);
  rowsPlan_c2r_ = fftw_plan_many_dft_c2r(xRank, xN, xHowmany, rowsBufC_, xOnembed, xOstride, xOdist,
    rowsBufR_, xInembed, xIstride, xIdist, FFTW_ESTIMATE);

  // Columns setup
  int yRank = 1;
  int yN[] = {static_cast<int>(ny_)};
  int yHowmany = nkPerTask_[myrank_]*nvz_*2;
  int *yInembed = NULL;
  const int yIstride = 1;
  const int yIdist = static_cast<int>(ny_);
  int *yOnembed = NULL;
  const int yOstride = 1;
  const int yOdist = static_cast<int>(nl_);
  colsBufR_ = fftw_alloc_real(yHowmany*ny_);
  colsBufC_ = fftw_alloc_complex(yHowmany*nl_);
  colsPlan_r2c_ = fftw_plan_many_dft_r2c(yRank, yN, yHowmany, colsBufR_, yInembed, yIstride, yIdist,
    colsBufC_, yOnembed, yOstride, yOdist, FFTW_ESTIMATE);
  colsPlan_c2r_ = fftw_plan_many_dft_c2r(yRank, yN, yHowmany, colsBufC_, yOnembed, yOstride, yOdist,
    colsBufR_, yInembed, yIstride, yIdist, FFTW_ESTIMATE);

  oops::Log::trace() << classname() << "::setupFFT done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransform::addSpectralCoefficient(const size_t & jk,
                                                const size_t & jl,
                                                const Quad & jq,
                                                const size_t & jsXDerivativeOffset,
                                                const size_t & jsYDerivativeOffset) {
  if (jk < ellips_.size()) {
    if (jl <= ellips_[jk]) {
      // Update spVec
      spElem e;
      e.jk = jk;
      e.jl = jl;
      e.jq = jq;
      e.jwGlb = static_cast<size_t>(kstar(jk, jl, nk_, nl_, nwGlb_)+jwGlbTol_);
      e.jsXDerivativeOffset = jsXDerivativeOffset;
      e.jsYDerivativeOffset = jsYDerivativeOffset;
      spVec_.push_back(e);

      // Check consistency between truncation and maximum total wavenumber
      ASSERT(e.jwGlb < nwGlb_);
    }
  }

  // Update spNorm
  ++spNormKL_[jk*nl_+jl];
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
