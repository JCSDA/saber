/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

class LayerRC : public LayerBase {
 public:
  static const std::string classname() {return "saber::fastlam::LayerRC";}

  // Constructor
  LayerRC(const FastLAMParametersBase & params,
          const oops::GeometryData & gdata,
          const std::string & myVar,
          const size_t & nx0,
          const size_t & ny0,
          const size_t & nz0) :
    LayerBase(params, gdata, myVar, nx0, ny0, nz0) {}
    ~LayerRC() = default;

  // Setups
  void setupParallelization() override;
  void setupNormalization() override;

  // Multiply square-root and adjoint
  size_t ctlVecSize() const override {return nxPerTask_[myrank_]*ny_*nz_;};
  void multiplySqrt(const atlas::Field &,
                    atlas::Field &,
                    const size_t &) const override;
  void multiplySqrtTrans(const atlas::Field &,
                         atlas::Field &,
                         const size_t &) const override;

 private:
  void print(std::ostream &) const override;

  // Transforms
  void redToRows(const atlas::Field &, atlas::Field &) const;
  void rowsToRed(const atlas::Field &, atlas::Field &) const;
  void rowsToCols(const atlas::Field &, atlas::Field &) const;
  void colsToRows(const atlas::Field &, atlas::Field &) const;

  // Convolutions
  void rowsConvolutionTL(atlas::Field &) const;
  void rowsConvolutionAD(atlas::Field &) const;
  void colsConvolutionTL(atlas::Field &) const;
  void colsConvolutionAD(atlas::Field &) const;
  void vertConvolution(atlas::Field &) const;

  // Normalizations
  void rowsNormalization(atlas::Field &) const;
  void colsNormalization(atlas::Field &) const;
  void vertNormalization(atlas::Field &) const;

  // Multiply square-root on reduced grid
  void multiplyRedSqrt(const atlas::Field &, atlas::Field &) const;
  void multiplyRedSqrtTrans(const atlas::Field &, atlas::Field &) const;

  // Sizes
  std::vector<size_t> nxPerTask_;
  std::vector<size_t> nyPerTask_;
  std::vector<size_t> nxStart_;
  std::vector<size_t> nyStart_;
  std::vector<size_t> nxEnd_;
  std::vector<size_t> nyEnd_;

  // Rows <=> reduced grid
  size_t xSendSize_;
  size_t rRecvSize_;
  std::vector<int> xTask_;
  std::vector<int> xOffset_;
  std::vector<int> xIndex_i_;
  std::vector<int> xIndex_j_;
  std::vector<int> xSendCounts_;
  std::vector<int> xSendDispls_;
  std::vector<int> rRecvCounts_;
  std::vector<int> rRecvDispls_;

  // Columns <=> rows
  size_t ySendSize_;
  size_t xRecvSize_;
  std::vector<std::vector<int>> yTask_;
  std::vector<std::vector<int>> yOffset_;
  std::vector<int> ySendCounts_;
  std::vector<int> ySendDispls_;
  std::vector<int> xRecvCounts_;
  std::vector<int> xRecvDispls_;
  std::vector<int> yIndex_i_;
  std::vector<int> yIndex_j_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
