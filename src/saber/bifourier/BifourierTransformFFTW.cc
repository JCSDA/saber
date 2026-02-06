/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierTransformFFTW.h"

#include <utility>

#include "saber/bifourier/BifourierUtilities.h"

using atlas::array::make_datatype;
using atlas::array::make_indexview;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static BifourierTransformMaker<BifourierTransformFFTW> makerFFTW_("fftw");

// -----------------------------------------------------------------------------

BifourierTransformFFTW::BifourierTransformFFTW(const oops::GeometryData & gdata,
                                               const std::string & gridUid,
                                               const oops::Variables & activeVars,
                                               const BifourierTransformParameters & params) :
  BifourierTransformBase(gdata, gridUid, activeVars, params) {
  oops::Log::trace() << classname() << "::BifourierTransformFFTW starting" << std::endl;

  // Setup global spectral space parameters
  setupGlobalSpectralSpace();

  // Setup parallelization, initial step [base method]
  setupParallelizationInit();

  // Backend-specific setup
  setupBackend();

  // Setup parallelization, final step [base method]
  setupParallelizationFinal();

  // Setup local spectral space
  setupLocalSpectralSpace();

  oops::Log::trace() << classname() << "::BifourierTransformFFTW done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierTransformFFTW::~BifourierTransformFFTW() {
  oops::Log::trace() << classname() << "::~BifourierTransformFFTW starting" << std::endl;

  // Delete FFTW-related data
  fftw_destroy_plan(rowsPlan_r2c_);
  fftw_destroy_plan(rowsPlan_c2r_);
  fftw_destroy_plan(colsPlan_r2c_);
  fftw_destroy_plan(colsPlan_c2r_);
  fftw_free(rowsBufR_);
  fftw_free(rowsBufC_);
  fftw_free(colsBufR_);
  fftw_free(colsBufC_);

  oops::Log::trace() << classname() << "::~BifourierTransformFFTW done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformFFTW::gp2sp(std::vector<double> & recvVec,
                                   atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::gp2sp starting" << std::endl;

  // Create send vector
  std::vector<double> sendVec;

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

  oops::Log::trace() << classname() << "::gp2sp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformFFTW::sp2gp(std::vector<double> & recvVec,
                                   const atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::sp2gp starting" << std::endl;

  // Create send vector
  std::vector<double> sendVec;

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

  oops::Log::trace() << classname() << "::sp2gp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformFFTW::setupBackend() {
  oops::Log::trace() << classname() << "::setupBackend starting" << std::endl;

  // Split in y direction
  nyPerTask_.resize(comm_.size(), 0);
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
  nkPerTask_.resize(comm_.size(), 0);
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

  // Rows <=> grid

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Index fields views
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  const auto indexIView = make_indexview<int, 1>(fs.index_i());
  const auto indexJView = make_indexview<int, 1>(fs.index_j());

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
        if (static_cast<size_t>(indexJView(jnode)) >= nyStart[jt] &&
          static_cast<size_t>(indexJView(jnode)) <= nyEnd[jt]) {
          rowsSendTask[jgr] = jt;
          rowsSendOffset[jgr] = rowsSendOffsetPerTask[jt];
          ++rowsSendOffsetPerTask[jt];
        }
      }
      ++jgr;
    }
  }

  // RecvCounts
  gridRecvCounts_.resize(comm_.size(), 0);
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
      gridRecvIndex_x[jgr] = indexIView(jnode);
      gridRecvIndex_y[jgr] = indexJView(jnode);
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
  rowsRecvCounts_.resize(comm_.size(), 0);
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

  // Equal chunks <=> columns

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
  colsRecvCounts_.resize(comm_.size(), 0);
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
  rowsPlan_r2c_ = fftw_plan_many_dft_r2c(xRank, xN, xHowmany, rowsBufR_, xInembed, xIstride,
    xIdist, rowsBufC_, xOnembed, xOstride, xOdist, FFTW_ESTIMATE);
  rowsPlan_c2r_ = fftw_plan_many_dft_c2r(xRank, xN, xHowmany, rowsBufC_, xOnembed, xOstride,
    xOdist, rowsBufR_, xInembed, xIstride, xIdist, FFTW_ESTIMATE);

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
  colsPlan_r2c_ = fftw_plan_many_dft_r2c(yRank, yN, yHowmany, colsBufR_, yInembed, yIstride,
    yIdist, colsBufC_, yOnembed, yOstride, yOdist, FFTW_ESTIMATE);
  colsPlan_c2r_ = fftw_plan_many_dft_c2r(yRank, yN, yHowmany, colsBufC_, yOnembed, yOstride,
    yOdist, colsBufR_, yInembed, yIstride, yIdist, FFTW_ESTIMATE);

  oops::Log::trace() << classname() << "::setupBackend done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
