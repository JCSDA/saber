/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <string>
#include <vector>

#include "ectrans/transi.h"

#include "saber/bifourier/BifourierTransformBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransformECTRANS : public BifourierTransformBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransformECTRANS";}

  // Constructor
  BifourierTransformECTRANS(const oops::GeometryData & gdata,
                            const std::string & gridUid,
                            const oops::Variables & activeVars,
                            const BifourierTransformParameters & params);

  // Destructor
  ~BifourierTransformECTRANS();

  // Forward FFT
  void gp2sp(std::vector<double> &,
             atlas::util::Metadata &) const;

  // Inverse FFT
  void sp2gp(std::vector<double> &,
             const atlas::util::Metadata &) const;

  // Forward FFT adjoint
  void gp2spAdj(std::vector<double> &,
                const atlas::util::Metadata &) const;

  // Inverse FFT adjoint
  void sp2gpAdj(std::vector<double> &,
                atlas::util::Metadata &) const;

  // Communication vectors
  const size_t & recvSize() const
    {return transSpRecvSize_;}
  const std::vector<int> & recvCounts() const
    {return transSpRecvCounts_;}
  const std::vector<int> & recvDispls() const
    {return transSpRecvDispls_;}

 private:
  // ECTRANS backend

  // Conversion from Fortran to C arrays
  const int f2c = -1;

  // Transform structure
  mutable struct Trans_t trans_;

  // Trans <=> grid
  size_t transGpSize_;
  size_t transGpSendSize_;
  std::vector<int> transGpSendIndex_;
  std::vector<int> transGpSendCounts_;
  std::vector<int> transGpSendDispls_;

  // Equal chunks <=> trans
  size_t transSpSize_;
  size_t transSpRecvSize_;
  std::vector<int> transSpRecvIndex_;
  std::vector<std::vector<int>> transSpJq_;
  std::vector<int> transSpRecvCounts_;
  std::vector<int> transSpRecvDispls_;

  // Private methods

  // Setup transform
  void setupBackend();
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
