/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fftw3.h>

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"

#include "saber/fastlam/FastLAMParametersBase.h"
#include "saber/fastlam/LayerBase.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class LayerSpec : public LayerBase {
 public:
  static const std::string classname() {return "saber::fastlam::LayerSpec";}

  // Constructor
  LayerSpec(const FastLAMParametersBase & params,
            const eckit::LocalConfiguration & fieldsMetaData,
            const oops::GeometryData & gdata,
            const std::string & myGroup,
            const std::vector<std::string> & myVars,
            const size_t & nx0,
            const size_t & ny0,
            const size_t & nz0) :
    LayerBase(params, fieldsMetaData, gdata, myGroup, myVars, nx0, ny0, nz0) {}
    ~LayerSpec() = default;

  // Setups
  void setupParallelization() override;
  void extractConvolution(const size_t &,
                          const size_t &,
                          std::vector<double> &,
                          std::vector<double> &) override;

  // Multiply square-root and adjoint
  size_t ctlVecSize() const override {return nxPerTask_[myrank_]*nyExt_*nz_;};
  void multiplySqrt(const atlas::Field &,
                    atlas::Field &,
                    const size_t &) const override;
  void multiplySqrtTrans(const atlas::Field &,
                         atlas::Field &,
                         const size_t &) const override;

 private:
  void print(std::ostream &) const override;

  // Transforms (including extension)
  void redToRows(const atlas::Field &, atlas::Field &) const;
  void rowsToRed(const atlas::Field &, atlas::Field &) const;
  void rowsToCols(const atlas::Field &, atlas::Field &) const;
  void colsToRows(const atlas::Field &, atlas::Field &) const;

  // Convolutions
  void rowsConvolution(atlas::Field &) const;
  void colsConvolution(atlas::Field &) const;
  void vertConvolution(atlas::Field &) const;

  // Normalizations
  void rowsNormalization(atlas::Field &) const;
  void colsNormalization(atlas::Field &) const;
  void vertNormalization(atlas::Field &) const;

  // Multiply square-root on reduced grid
  void multiplyRedSqrt(const atlas::Field &, atlas::Field &) const;
  void multiplyRedSqrtTrans(const atlas::Field &, atlas::Field &) const;

  // Sizes
  size_t nxExt_;
  size_t nyExt_;
  std::vector<size_t> nxPerTask_;
  std::vector<size_t> nyPerTask_;
  std::vector<size_t> nxStart_;
  std::vector<size_t> nyStart_;
  std::vector<size_t> nxEnd_;
  std::vector<size_t> nyEnd_;

  // Rows <=> reduced grid
  size_t xSendSize_;
  size_t rRecvSize_;
  std::vector<int> xSendTask_;
  std::vector<int> xSendOffset_;
  std::vector<int> xSendIndex_x_;
  std::vector<int> xSendIndex_y_;
  std::vector<int> xSendCounts_;
  std::vector<int> xSendDispls_;
  std::vector<int> rRecvCounts_;
  std::vector<int> rRecvDispls_;

  // Columns <=> rows
  size_t ySendSize_;
  size_t xRecvSize_;
  std::vector<int> ySendTask_;
  std::vector<int> ySendOffset_;
  std::vector<int> ySendCounts_;
  std::vector<int> ySendDispls_;
  std::vector<int> xRecvCounts_;
  std::vector<int> xRecvDispls_;
  std::vector<int> ySendIndex_x_;
  std::vector<int> ySendIndex_y_;

  // Rows FFT
  size_t nk_;
  fftw_plan xPlan_r2c_;
  fftw_plan xPlan_c2r_;
  double *xBufR_;
  fftw_complex *xBufC_;
  double xNormFFT_;
  std::vector<double> xSpecStdDev_;

  // Columns FFT
  size_t nl_;
  fftw_plan yPlan_r2c_;
  fftw_plan yPlan_c2r_;
  double *yBufR_;
  fftw_complex *yBufC_;
  double yNormFFT_;
  std::vector<double> ySpecStdDev_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
