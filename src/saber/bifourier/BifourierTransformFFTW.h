/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <fftw3.h>

#include <string>
#include <vector>

#include "saber/bifourier/BifourierTransformBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransformFFTW : public BifourierTransformBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransformFFTW";}

  // Constructor
  BifourierTransformFFTW(const oops::GeometryData &,
                         const std::string &,
                         const oops::Variables &,
                         const BifourierTransformParameters &);

  // Destructor
  ~BifourierTransformFFTW();

  // Forward FFT
  void gp2sp(std::vector<double> &,
             atlas::util::Metadata &) const;

  // Inverse FFT
  void sp2gp(std::vector<double> &,
             const atlas::util::Metadata &) const;

  // Forward FFT adjoint
  void gp2spAdj(std::vector<double> & recvVec,
                const atlas::util::Metadata & metadata) const
    {sp2gp(recvVec, metadata);}

  // Inverse FFT adjoint
  void sp2gpAdj(std::vector<double> & recvVec,
                atlas::util::Metadata & metadata) const
    {gp2sp(recvVec, metadata);}

  // Forward FFT normalization
  double gp2spNorm(const size_t & js) const
    {return normFFT_;}

  // Inverse FFT normalization
  double sp2gpNorm(const size_t & js) const
    {return 1.0;}

  // Forward FFT adjoint normalization
  double gp2spAdjNorm(const size_t & js) const
    {return normFFT_/spNorm(js);}

  // Inverse FFT adjoint normalization
  double sp2gpAdjNorm(const size_t & js) const
    {return spNorm(js);}

  // Communication vectors
  const size_t & recvSize() const
    {return colsRecvSize_;}
  const std::vector<int> & recvCounts() const
    {return colsRecvCounts_;}
  const std::vector<int> & recvDispls() const
    {return colsRecvDispls_;}

 private:
  // Sizes
  std::vector<size_t> nyPerTask_;
  std::vector<size_t> nkPerTask_;

  // Rows <=> grid
  size_t rowsSendSize_;
  std::vector<int> rowsSendIndex_;
  std::vector<int> rowsSendCounts_;
  std::vector<int> rowsSendDispls_;

  // Columns <=> rows
  size_t colsSendSize_;
  size_t rowsRecvSize_;
  std::vector<int> rowsRecvIndex_;
  std::vector<int> colsSendIndex_;
  std::vector<int> colsSendCounts_;
  std::vector<int> colsSendDispls_;
  std::vector<int> rowsRecvCounts_;
  std::vector<int> rowsRecvDispls_;

  // Equal chunks <=> columns
  size_t colsRecvSize_;
  std::vector<int> colsRecvIndex_;
  std::vector<std::vector<int>> colsJq_;
  std::vector<int> colsRecvCounts_;
  std::vector<int> colsRecvDispls_;

  // Rows FFT
  fftw_plan rowsPlan_r2c_;
  fftw_plan rowsPlan_c2r_;
  double *rowsBufR_ = nullptr;
  fftw_complex *rowsBufC_ = nullptr;

  // Columns FFT
  fftw_plan colsPlan_r2c_;
  fftw_plan colsPlan_c2r_;
  double *colsBufR_ = nullptr;
  fftw_complex *colsBufC_ = nullptr;

  // Private methods

  // Setup backend
  void setupBackend();
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
