/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <tuple>
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

class LayerHalo : public LayerBase {
 public:
  static const std::string classname() {return "saber::fastlam::LayerHalo";}

  // Constructor
  LayerHalo(const FastLAMParametersBase & params,
            const oops::GeometryData & gdata,
            const std::string & myVar,
            const size_t & nx0,
            const size_t & ny0,
            const size_t & nz0) :
    LayerBase(params, gdata, myVar, nx0, ny0, nz0) {}
  ~LayerHalo() = default;

  // Setups
  void setupParallelization() override;
  void setupNormalization() override;

  // Multiply square-root and adjoint
  size_t ctlVecSize() const override {return rSize_*nz_;};
  void multiplySqrt(const atlas::Field &,
                    atlas::Field &,
                    const size_t &) const override;
  void multiplySqrtTrans(const atlas::Field &,
                         atlas::Field &,
                         const size_t &) const override;

 private:
  void print(std::ostream &) const override;

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
  void multiplyRedSqrt(atlas::Field &) const;
  void multiplyRedSqrtTrans(atlas::Field &) const;

  // Rows <=> reduced grid
  size_t xcSize_;
  std::vector<std::tuple<int, int, double>> xcOperations_;
  size_t xcRecvSize_;
  std::vector<int> xcRecvCounts_;
  std::vector<int> xcRecvDispls_;
  size_t xcSendSize_;
  std::vector<int> xcSendCounts_;
  std::vector<int> xcSendDispls_;
  std::vector<int> xcSendMapping_;

  // Columns <=> rows
  size_t ycSize_;
  std::vector<std::tuple<int, int, double>> ycOperations_;
  size_t ycRecvSize_;
  std::vector<int> ycRecvCounts_;
  std::vector<int> ycRecvDispls_;
  size_t ycSendSize_;
  std::vector<int> ycSendCounts_;
  std::vector<int> ycSendDispls_;
  std::vector<int> ycSendMapping_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
