/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierTransformECTRANS.h"

#include <algorithm>
#include <limits>
#include <utility>

#include "saber/bifourier/BifourierUtilities.h"

using atlas::array::make_datatype;
using atlas::array::make_indexview;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static BifourierTransformMaker<BifourierTransformECTRANS> makerECTRANS_("ectrans");

// -----------------------------------------------------------------------------

BifourierTransformECTRANS::BifourierTransformECTRANS(const oops::GeometryData & gdata,
                                                     const std::string & gridUid,
                                                     const oops::Variables & activeVars,
                                                     const BifourierTransformParameters & params) :
  BifourierTransformBase(gdata, gridUid, activeVars, params) {
  oops::Log::trace() << classname() << "::BifourierTransformECTRANS starting" << std::endl;

  // Setup global spectral space parameters
  setupGlobalSpectralSpace();

  // Setup parallelization, initial step [base method]
  setupParallelizationInit();

  // Setup parallelization, backend-specific
  setupBackend();

  // Setup parallelization, final step [base method]
  setupParallelizationFinal();

  // Setup local spectral space
  setupLocalSpectralSpace();

  oops::Log::trace() << classname() << "::BifourierTransformECTRANS done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierTransformECTRANS::~BifourierTransformECTRANS() {
  oops::Log::trace() << classname() << "::~BifourierTransformECTRANS starting" << std::endl;

  // Delete and finalize transform structure
  trans_delete(&trans_);
  trans_finalize();

  oops::Log::trace() << classname() << "::~BifourierTransformECTRANS done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformECTRANS::gp2sp(std::vector<double> & recvVec,
                                      atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::gp2sp starting" << std::endl;

  // Check metadata for wind transform
  const size_t nvordiv = metadata.getInt("nvordiv", 0);

  // Number of scalar fields
  const size_t nscalar = nvz_-2*nvordiv;

  // Allocate buffers
  std::vector<double> rgp(transGpSize_*nvz_);
  std::vector<double> rspScalar;
  std::vector<double> rspVor;
  std::vector<double> rspDiv;
  std::vector<double> rMeanU;
  std::vector<double> rMeanV;
  rspScalar.resize(transSpSize_*nscalar);
  rspVor.resize(transSpSize_*nvordiv);
  rspDiv.resize(transSpSize_*nvordiv);
  rMeanU.resize(nvordiv);
  rMeanV.resize(nvordiv);

  // Setup direct transform
  struct DirTrans_t dirtrans;
  dirtrans = new_dirtrans(&trans_);
  dirtrans.nscalar = nscalar;
  dirtrans.nvordiv = nvordiv;
  dirtrans.rspscalar = rspScalar.data();
  dirtrans.rspvor = rspVor.data();
  dirtrans.rspdiv = rspDiv.data();
  dirtrans.rmeanu = rMeanU.data();
  dirtrans.rmeanv = rMeanV.data();
  dirtrans.rgp = rgp.data();

  // Communication
  std::vector<double> sendVec(transGpSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    sendVec.data(), transGpSendCounts_.data(), transGpSendDispls_.data());

  // Reserialize
  for (size_t jrs = 0; jrs < transGpSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      const size_t jf = jvz*transGpSendSize_ + transGpSendIndex_[jrs];

      // Copy data
      rgp[jf] = sendVec[jrsv];
    }
  }

  // Compute direct transform
  trans_dirtrans(&dirtrans);

  if (nvordiv > 0) {
    // Save mean wind
    const double windNorm = std::sqrt(normFFT_);
    for (size_t jz = 0; jz < nvordiv; ++jz) {
      rMeanU[jz] *= windNorm;
      rMeanV[jz] *= windNorm;
    }
    metadata.set("uMeanProfile", rMeanU);
    metadata.set("vMeanProfile", rMeanV);
  }

  // Reserialize
  recvVec.resize(transSpRecvSize_*nvz_);
  size_t jj = 0;
  size_t jfBase = 0;
  for (int jk = 0; jk < trans_.nump; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
      for (size_t jc = 0; jc < transSpJq_[jk*nl_+jl].size(); ++jc) {
        // Get jq
        const size_t jq = transSpJq_[jk*nl_+jl][jc];

        // Initialize variable and level counter
        size_t jvz = 0;

        // Vorticity
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          recvVec[jcrv] = rspVor[jf];

          // Update counter
          ++jvz;
        }

        // Divergence
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          recvVec[jcrv] = rspDiv[jf];

          // Update counter
          ++jvz;
        }

        // Scalars
        for (size_t jz = 0; jz < nscalar; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nscalar + jq*nscalar + jz;

          // Copy data
          recvVec[jcrv] = rspScalar[jf];

          // Update counter
          ++jvz;
        }

        // Check counter
        ASSERT(jvz == nvz_);

        // Update communication vector index
        ++jj;
      }

      // Update FFT vector base index
      if (transSpJq_[jk*nl_+jl].size() > 0) {
        ++jfBase;
      }
    }
  }

  oops::Log::trace() << classname() << "::gp2sp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformECTRANS::sp2gp(std::vector<double> & recvVec,
                                      const atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::sp2gp starting" << std::endl;

  // Check metadata for wind transform
  const size_t nvordiv = metadata.getInt("nvordiv", 0);

  // Number of scalar fields
  const size_t nscalar = nvz_-2*nvordiv;

  // Allocate buffers
  std::vector<double> rgp(transGpSize_*nvz_);
  std::vector<double> rspScalar;
  std::vector<double> rspVor;
  std::vector<double> rspDiv;
  std::vector<double> rMeanU;
  std::vector<double> rMeanV;
  rspScalar.resize(transSpSize_*nscalar, 0.0);
  rspVor.resize(transSpSize_*nvordiv, 0.0);
  rspDiv.resize(transSpSize_*nvordiv, 0.0);
  rMeanU.resize(nvordiv);
  rMeanV.resize(nvordiv);

  // Setup inverse transform
  struct InvTrans_t invtrans;
  invtrans = new_invtrans(&trans_);
  invtrans.nscalar = nscalar;
  invtrans.nvordiv = nvordiv;
  invtrans.rspscalar = rspScalar.data();
  invtrans.rspvor = rspVor.data();
  invtrans.rspdiv = rspDiv.data();
  invtrans.rmeanu = rMeanU.data();
  invtrans.rmeanv = rMeanV.data();
  invtrans.rgp = rgp.data();

  // Reserialize
  size_t jj = 0;
  size_t jfBase = 0;
  for (int jk = 0; jk < trans_.nump; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
      for (size_t jc = 0; jc < transSpJq_[jk*nl_+jl].size(); ++jc) {
        // Get jq
        const size_t jq = transSpJq_[jk*nl_+jl][jc];

        // Initialize variable and level counter
        size_t jvz = 0;

        // Vorticity
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          rspVor[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Divergence
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          rspDiv[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Scalars
        for (size_t jz = 0; jz < nscalar; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nscalar + jq*nscalar + jz;

          // Copy data
          rspScalar[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Check counter
        ASSERT(jvz == nvz_);

        // Update communication vector index
        ++jj;
      }

      // Update FFT vector base index
      if (transSpJq_[jk*nl_+jl].size() > 0) {
        ++jfBase;
      }
    }
  }

  if (metadata.has("uMeanProfile") && metadata.has("vMeanProfile")) {
    // Copy mean wind
    const std::vector<double> uMeanProfile = metadata.getDoubleVector("uMeanProfile");
    const std::vector<double> vMeanProfile = metadata.getDoubleVector("vMeanProfile");
    for (size_t jz = 0; jz < nvordiv; ++jz) {
      rMeanU[jz] = uMeanProfile[jz];
      rMeanV[jz] = vMeanProfile[jz];
    }
  } else {
    // Set mean wind to zero
    std::fill(rMeanU.begin(), rMeanU.end(), 0.0);
    std::fill(rMeanV.begin(), rMeanV.end(), 0.0);
  }

  // Compute inverse transform
  trans_invtrans(&invtrans);

  // Reserialize
  std::vector<double> sendVec(transGpSendSize_*nvz_);
  for (size_t jrs = 0; jrs < transGpSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      const size_t jf = jvz*transGpSendSize_ + transGpSendIndex_[jrs];

      // Copy data
      sendVec[jrsv] = rgp[jf];
    }
  }

  // Communication
  recvVec.resize(gridRecvSize_*nvz_);
  comm_.allToAllv(sendVec.data(), transGpSendCounts_.data(), transGpSendDispls_.data(),
    recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data());

  oops::Log::trace() << classname() << "::sp2gp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformECTRANS::gp2spAdj(std::vector<double> & recvVec,
                                         const atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::gp2spAdj starting" << std::endl;

  // Check metadata for wind transform
  const size_t nvordiv = metadata.getInt("nvordiv", 0);

  // Number of scalar fields
  const size_t nscalar = nvz_-2*nvordiv;

  // Allocate buffers
  std::vector<double> rgp(transGpSize_*nvz_);
  std::vector<double> rspScalar;
  std::vector<double> rspVor;
  std::vector<double> rspDiv;
  std::vector<double> rMeanU;
  std::vector<double> rMeanV;
  rspScalar.resize(transSpSize_*nscalar, 0.0);
  rspVor.resize(transSpSize_*nvordiv, 0.0);
  rspDiv.resize(transSpSize_*nvordiv, 0.0);
  rMeanU.resize(nvordiv);
  rMeanV.resize(nvordiv);

  // Setup inverse transform
  struct DirTransAdj_t dirtrans_adj;
  dirtrans_adj = new_dirtrans_adj(&trans_);
  dirtrans_adj.nscalar = nscalar;
  dirtrans_adj.nvordiv = nvordiv;
  dirtrans_adj.rspscalar = rspScalar.data();
  dirtrans_adj.rspvor = rspVor.data();
  dirtrans_adj.rspdiv = rspDiv.data();
  dirtrans_adj.rmeanu = rMeanU.data();
  dirtrans_adj.rmeanv = rMeanV.data();
  dirtrans_adj.rgp = rgp.data();

  // Reserialize
  size_t jj = 0;
  size_t jfBase = 0;
  for (int jk = 0; jk < trans_.nump; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
      for (size_t jc = 0; jc < transSpJq_[jk*nl_+jl].size(); ++jc) {
        // Get jq
        const size_t jq = transSpJq_[jk*nl_+jl][jc];

        // Initialize variable and level counter
        size_t jvz = 0;

        // Vorticity
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

           // Copy data
          rspVor[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Divergence
        for (size_t jz = 0; jz < nvordiv; ++jz) {
           // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          rspDiv[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Scalars
        for (size_t jz = 0; jz < nscalar; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nscalar + jq*nscalar + jz;

          // Copy data
          rspScalar[jf] = recvVec[jcrv];

          // Update counter
          ++jvz;
        }

        // Check counter
        ASSERT(jvz == nvz_);

        // Update communication vector index
        ++jj;
      }

      // Update FFT vector base index
      if (transSpJq_[jk*nl_+jl].size() > 0) {
        ++jfBase;
      }
    }
  }

  if (metadata.has("uMeanProfile") && metadata.has("vMeanProfile")) {
    // Copy mean wind
    const std::vector<double> uMeanProfile = metadata.getDoubleVector("uMeanProfile");
    const std::vector<double> vMeanProfile = metadata.getDoubleVector("vMeanProfile");
    for (size_t jz = 0; jz < nvordiv; ++jz) {
      rMeanU[jz] = uMeanProfile[jz];
      rMeanV[jz] = vMeanProfile[jz];
    }
  } else {
    // Set mean wind to zero
    std::fill(rMeanU.begin(), rMeanU.end(), 0.0);
    std::fill(rMeanV.begin(), rMeanV.end(), 0.0);
  }

  // Compute direct adjoint transform
  trans_dirtrans_adj(&dirtrans_adj);

  // Reserialize
  std::vector<double> sendVec(transGpSendSize_*nvz_);
  for (size_t jrs = 0; jrs < transGpSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      const size_t jf = jvz*transGpSendSize_ + transGpSendIndex_[jrs];

      // Copy data
      sendVec[jrsv] = rgp[jf];
    }
  }

  // Communication
  recvVec.resize(gridRecvSize_*nvz_);
  comm_.allToAllv(sendVec.data(), transGpSendCounts_.data(), transGpSendDispls_.data(),
    recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data());

  oops::Log::trace() << classname() << "::gp2spAdj done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformECTRANS::sp2gpAdj(std::vector<double> & recvVec,
                                         atlas::util::Metadata & metadata) const {
  oops::Log::trace() << classname() << "::sp2gpAdj starting" << std::endl;

  // Check metadata for wind transform
  const size_t nvordiv = metadata.getInt("nvordiv", 0);

  // Number of scalar fields
  const size_t nscalar = nvz_-2*nvordiv;

  // Allocate buffers
  std::vector<double> rgp(transGpSize_*nvz_);
  std::vector<double> rspScalar;
  std::vector<double> rspVor;
  std::vector<double> rspDiv;
  std::vector<double> rMeanU;
  std::vector<double> rMeanV;
  rspScalar.resize(transSpSize_*nscalar);
  rspVor.resize(transSpSize_*nvordiv);
  rspDiv.resize(transSpSize_*nvordiv);
  rMeanU.resize(nvordiv);
  rMeanV.resize(nvordiv);

  // Setup direct transform
  struct InvTransAdj_t invtrans_adj;
  invtrans_adj = new_invtrans_adj(&trans_);
  invtrans_adj.nscalar = nscalar;
  invtrans_adj.nvordiv = nvordiv;
  invtrans_adj.rspscalar = rspScalar.data();
  invtrans_adj.rspvor = rspVor.data();
  invtrans_adj.rspdiv = rspDiv.data();
  invtrans_adj.rmeanu = rMeanU.data();
  invtrans_adj.rmeanv = rMeanV.data();
  invtrans_adj.rgp = rgp.data();

  // Communication
  std::vector<double> sendVec(transGpSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    sendVec.data(), transGpSendCounts_.data(), transGpSendDispls_.data());

  // Reserialize
  for (size_t jrs = 0; jrs < transGpSendSize_; ++jrs) {
    for (size_t jvz = 0; jvz < nvz_; ++jvz) {
      // Communication vector index
      const size_t jrsv = jrs*nvz_ + jvz;

      // FFT vector index
      const size_t jf = jvz*transGpSendSize_ + transGpSendIndex_[jrs];

      // Copy data
      rgp[jf] = sendVec[jrsv];
    }
  }

  // Compute inverse adjoint transform
  trans_invtrans_adj(&invtrans_adj);

  if (nvordiv > 0) {
    // Save mean wind
    const double windNorm = std::sqrt(normFFT_);
    for (size_t jz = 0; jz < nvordiv; ++jz) {
      rMeanU[jz] *= windNorm;
      rMeanV[jz] *= windNorm;
    }
    metadata.set("uMeanProfile", rMeanU);
    metadata.set("vMeanProfile", rMeanV);
  }

  // Reserialize
  recvVec.resize(transSpRecvSize_*nvz_);
  size_t jj = 0;
  size_t jfBase = 0;
  for (int jk = 0; jk < trans_.nump; ++jk) {
    for (size_t jl = 0; jl < nl_; ++jl) {
      for (size_t jc = 0; jc < transSpJq_[jk*nl_+jl].size(); ++jc) {
        // Get jq
        const size_t jq = transSpJq_[jk*nl_+jl][jc];

        // Initialize variable and level counter
        size_t jvz = 0;

        // Vorticity
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          recvVec[jcrv] = rspVor[jf];

          // Update counter
          ++jvz;
        }

        // Divergence
        for (size_t jz = 0; jz < nvordiv; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nvordiv + jq*nvordiv + jz;

          // Copy data
          recvVec[jcrv] = rspDiv[jf];

          // Update counter
          ++jvz;
        }

        // Scalars
        for (size_t jz = 0; jz < nscalar; ++jz) {
          // Communication vector index
          const size_t jcrv = transSpRecvIndex_[jj] + jvz;

          // FFT vector index
          const size_t jf = jfBase*4*nscalar + jq*nscalar + jz;

          // Copy data
          recvVec[jcrv] = rspScalar[jf];

          // Update counter
          ++jvz;
        }

        // Check counter
        ASSERT(jvz == nvz_);

        // Update communication vector index
        ++jj;
      }

      // Update FFT vector base index
      if (transSpJq_[jk*nl_+jl].size() > 0) {
        ++jfBase;
      }
    }
  }

  oops::Log::trace() << classname() << "::sp2gpAdj done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformECTRANS::setupBackend() {
  oops::Log::trace() << classname() << "::setupBackend starting" << std::endl;

  // Configure transforms
  trans_use_mpi(true);
  trans_set_leq_regions(false);
  const int nprgpew = std::min(1, static_cast<int>(std::sqrt(static_cast<double>(comm_.size()))));
  trans_set_nprgpew(nprgpew);

  // Setup transform structure
  trans_new(&trans_);
  trans_set_resol_lam(&trans_, nx_, ny_, dx_, dy_);
  trans_set_trunc_lam(&trans_, M_, N_);
  trans_setup(&trans_);
  trans_inquire(&trans_, "nlstlat,nfrstlat,nsta,nonl,nmyms");

  // Copy transform structure sizes
  transGpSize_ = trans_.ngptot;
  transSpSize_ = trans_.nspec2;

  // Trans <=> grid

  // Define x/y indices from trans grid-point distribution (see IFS code for details...)
  std::vector<size_t> index_x(transGpSize_);
  std::vector<size_t> index_y(transGpSize_);
  size_t xMin = std::numeric_limits<size_t>::max();
  size_t xMax = std::numeric_limits<size_t>::min();
  size_t yMin = std::numeric_limits<size_t>::max();
  size_t yMax = std::numeric_limits<size_t>::min();
  int irof = 0;
  const int istlat = 1;
  const int iendlat = trans_.nlstlat[trans_.my_region_NS+f2c]-trans_.nfrstloff;
  for (int jgl = istlat; jgl <= iendlat; ++jgl) {
    const int iglg = trans_.nfrstlat[trans_.my_region_NS+f2c]+jgl-istlat;
    const int nstaIndex = (trans_.nptrfloff+jgl-1)*trans_.n_regions_EW + trans_.my_region_EW;
    const int istlon = trans_.nsta[nstaIndex+f2c];
    const int iendlon = trans_.nsta[nstaIndex+f2c]+trans_.nonl[nstaIndex+f2c]-1;
    for (int jlon = istlon; jlon <= iendlon; ++jlon) {
      const int igp = trans_.nlon*(iglg-1)+jlon;
      index_y[irof] = (igp-1)/trans_.nlon+1;
      index_x[irof] = igp-(index_y[irof]-1)*trans_.nlon;
      xMin = std::min(xMin, index_x[irof]);
      xMax = std::max(xMax, index_x[irof]);
      yMin = std::min(yMin, index_y[irof]);
      yMax = std::max(yMax, index_y[irof]);
      ++irof;
    }
  }
  ASSERT(irof == static_cast<int>(transGpSize_));

  // Trans grid local index from x/y
  const size_t nx = xMax-xMin+1;
  const size_t ny = yMax-yMin+1;
  std::vector<std::vector<int>> gpIndex(nx);
  for (size_t jx = 0; jx < nx; ++jx) {
    gpIndex[jx].resize(ny, -1);
  }
  for (size_t jj = 0; jj < transGpSize_; ++jj) {
    gpIndex[index_x[jj]-xMin][index_y[jj]-yMin] = jj;
  }

  // Copy x/y/rank
  std::vector<int> nto(3, 1);
  std::vector<double> rgpg;
  if (comm_.rank() == 0) {
    rgpg.resize(3*trans_.ngptotg);
  }
  std::vector<double> rgp(3*transGpSize_);
  for (size_t jj = 0; jj < transGpSize_; ++jj) {
    rgp[0*transGpSize_+jj] = static_cast<double>(index_x[jj]);
    rgp[1*transGpSize_+jj] = static_cast<double>(index_y[jj]);
    rgp[2*transGpSize_+jj] = comm_.rank();
  }

  // Gather x/y/rank
  struct GathGrid_t gathgrid = new_gathgrid(&trans_);
  gathgrid.rgp  = rgp.data();
  gathgrid.rgpg = rgpg.data();
  gathgrid.nto  = nto.data();
  gathgrid.nfld = 3;
  trans_gathgrid(&gathgrid);

  // Get partition
  std::vector<int> partgp(trans_.ngptotg);
  if (comm_.rank() == 0) {
    ASSERT(trans_.ngptotg == trans_.ngptotg);
    for (int jj = 0; jj < trans_.ngptotg; ++jj) {
      const size_t ix = static_cast<size_t>(rgpg[0*trans_.ngptotg+jj])+f2c;
      const size_t iy = static_cast<size_t>(rgpg[1*trans_.ngptotg+jj])+f2c;
      const size_t it = static_cast<size_t>(rgpg[2*trans_.ngptotg+jj]);
      const size_t ixy = ix*ny_+iy;
      partgp[ixy] = it;
    }
  }

  // Broadcast partition
  comm_.broadcast(partgp, 0);

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
  std::vector<int> transGpSendTask(gridRecvSize_);
  std::vector<int> transGpSendOffset(gridRecvSize_);
  std::vector<int> transGpSendOffsetPerTask(comm_.size(), 0);
  size_t jgr = 0;
  for (size_t jnode = 0; jnode < nodes_; ++jnode) {
    if (ghostView(jnode) == 0) {
      const size_t ix = indexIView(jnode);
      const size_t iy = indexJView(jnode);
      const size_t ixy = ix*ny_+iy;
      const size_t it = partgp[ixy];
      transGpSendTask[jgr] = it;
      transGpSendOffset[jgr] = transGpSendOffsetPerTask[it];
      ++transGpSendOffsetPerTask[it];
      ++jgr;
    }
  }

  // RecvCounts
  gridRecvCounts_.resize(comm_.size(), 0);
  for (size_t jgr = 0; jgr < gridRecvSize_; ++jgr) {
    ++gridRecvCounts_[transGpSendTask[jgr]];
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
    transGpSendCounts_.push_back(rRecvCountsGlb_[jt*comm_.size()+myrank_]);
  }

  // Buffer size
  transGpSendSize_ = 0;
  for (const auto & n : transGpSendCounts_) transGpSendSize_ += n;

  // SendDispls
  transGpSendDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    transGpSendDispls_[jt] = static_cast<int>(jt ?
      transGpSendDispls_[jt-1] + transGpSendCounts_[jt-1] : 0);
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
  std::vector<int> transGpSendIndex_x(transGpSendSize_);
  std::vector<int> transGpSendIndex_y(transGpSendSize_);
  transGpSendIndex_.resize(transGpSendSize_);
  comm_.allToAllv(gridRecvIndex_x.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    transGpSendIndex_x.data(), transGpSendCounts_.data(), transGpSendDispls_.data());
  comm_.allToAllv(gridRecvIndex_y.data(), gridRecvCounts_.data(), gridRecvDispls_.data(),
    transGpSendIndex_y.data(), transGpSendCounts_.data(), transGpSendDispls_.data());
  for (size_t jrs = 0; jrs < transGpSendSize_; ++jrs) {
    transGpSendIndex_[jrs] = gpIndex[transGpSendIndex_x[jrs]-(xMin+f2c)]
      [transGpSendIndex_y[jrs]-(yMin+f2c)];
    ASSERT(transGpSendIndex_[jrs] >= 0);
    ASSERT(transGpSendIndex_[jrs] < static_cast<int>(transGpSize_));
  }

  // Effective index
  gridRecvIndex_.resize(gridRecvSize_);
  for (size_t jgr = 0; jgr < gridRecvSize_; ++jgr) {
    gridRecvIndex_[jgr] = (gridRecvDispls_[transGpSendTask[jgr]] + transGpSendOffset[jgr])*nvz_;
  }

  // Equal chunks <=> trans

  // Number of values on each destination task
  transSpRecvSize_ = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    for (int jump = 1; jump <= trans_.nump; ++jump) {
      if (static_cast<int>(jk) == trans_.nmyms[jump+f2c]) {
        ++transSpRecvSize_;
      }
    }
  }

  // Define destination task
  std::vector<int> eqchSendTask(transSpRecvSize_);
  std::vector<int> eqchSendOffset(transSpRecvSize_);
  std::vector<int> eqchSendOffsetPerTask(comm_.size(), 0);
  size_t jv = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    for (int jump = 1; jump <= trans_.nump; ++jump) {
      if (static_cast<int>(jk) == trans_.nmyms[jump+f2c]) {
        const size_t jt = spVec_[jsGlb].jt;
        eqchSendTask[jv] = jt;
        eqchSendOffset[jv] = eqchSendOffsetPerTask[jt];
        ++eqchSendOffsetPerTask[jt];
        ++jv;
      }
    }
  }

  // RecvCounts
  transSpRecvCounts_.resize(comm_.size(), 0);
  for (size_t jcr = 0; jcr < transSpRecvSize_; ++jcr) {
    ++transSpRecvCounts_[eqchSendTask[jcr]];
  }

  // RecvDispls
  transSpRecvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    transSpRecvDispls_[jt] = static_cast<int>(jt ?
      transSpRecvDispls_[jt-1] + transSpRecvCounts_[jt-1] : 0);
  }

  // Allgather RecvCounts
  eckit::mpi::Buffer<int> transSpRecvCountsBuffer(comm_.size());
  comm_.allGatherv(transSpRecvCounts_.begin(), transSpRecvCounts_.end(), transSpRecvCountsBuffer);
  std::vector<int> transSpRecvCountsGlb_ = std::move(transSpRecvCountsBuffer.buffer);

  // SendCounts
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    eqchSendCounts_.push_back(transSpRecvCountsGlb_[jt*comm_.size()+myrank_]);
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
  std::vector<size_t> transSpRecvOrder(transSpRecvSize_);
  std::iota(transSpRecvOrder.begin(), transSpRecvOrder.end(), 0);
  std::stable_sort(transSpRecvOrder.begin(), transSpRecvOrder.end(),
    [&](int i, int j){return eqchSendTask[i] < eqchSendTask[j];});
  std::vector<size_t> transSpRecvOrderInverse(transSpRecvSize_);
  for (size_t jcr = 0; jcr < transSpRecvSize_; ++jcr) {
    transSpRecvOrderInverse[transSpRecvOrder[jcr]] = jcr;
  }

  // Communicate indices
  std::vector<int> transSpRecvIndex(transSpRecvSize_);
  size_t jcr = 0;
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    for (int jump = 1; jump <= trans_.nump; ++jump) {
      if (static_cast<int>(jk) == trans_.nmyms[jump+f2c]) {
        transSpRecvIndex[transSpRecvOrderInverse[jcr]] = jsGlb;
        ++jcr;
      }
    }
  }
  eqchSendIndex_.resize(eqchSendSize_);
  comm_.allToAllv(transSpRecvIndex.data(), transSpRecvCounts_.data(), transSpRecvDispls_.data(),
    eqchSendIndex_.data(), eqchSendCounts_.data(), eqchSendDispls_.data());
  for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
    const auto it = std::find(sToSGlb_.begin(), sToSGlb_.end(), eqchSendIndex_[jes]);
    ASSERT(it != sToSGlb_.end());
    eqchSendIndex_[jes] = it-sToSGlb_.begin();
  }

  // Effective index
  transSpRecvIndex_.resize(transSpRecvSize_);
  for (size_t jcr = 0; jcr < transSpRecvSize_; ++jcr) {
    transSpRecvIndex_[jcr] = (transSpRecvDispls_[eqchSendTask[jcr]] + eqchSendOffset[jcr])*nvz_;
  }

  // Scale counts and displs for all levels
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    gridRecvCounts_[jt] *= nvz_;
    gridRecvDispls_[jt] *= nvz_;
    transGpSendCounts_[jt] *= nvz_;
    transGpSendDispls_[jt] *= nvz_;
    transSpRecvCounts_[jt] *= nvz_;
    transSpRecvDispls_[jt] *= nvz_;
    eqchSendCounts_[jt] *= nvz_;
    eqchSendDispls_[jt] *= nvz_;
  }

  // Number of component for each local jk,jl couple
  transSpJq_.resize(trans_.nump*nl_);
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    const size_t jk = spVec_[jsGlb].jk;
    for (int jump = 1; jump <= trans_.nump; ++jump) {
      if (static_cast<int>(jk) == trans_.nmyms[jump+f2c]) {
        const size_t jkPerTask = jump+f2c;
        const size_t jl = spVec_[jsGlb].jl;
        const size_t jq = spVec_[jsGlb].jq;
        transSpJq_[jkPerTask*nl_+jl].push_back(jq);
      }
    }
  }

  oops::Log::trace() << classname() << "::setupBackend done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
