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
#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"

#include "saber/fastlam/FastLAMParametersBase.h"
#include "saber/fastlam/InterpElement.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class Layer {
 public:
  static const std::string classname() {return "saber::fastlam::Layer";}

  typedef FastLAMParametersBase ParametersBase_;

  // Constructor
  Layer(const ParametersBase_ & params,
        const oops::GeometryData & gdata,
        const std::string & myVar,
        const size_t & nx0,
        const size_t & ny0,
        const size_t & nz0) :
    params_(params),
    gdata_(gdata),
    comm_(gdata_.comm()),
    myrank_(comm_.rank()),
    myVar_(myVar),
    nx0_(nx0),
    ny0_(ny0),
    nodes0_(gdata_.functionSpace().ghost().shape(0)),
    nz0_(nz0) {}

  // Setups
  void setupVerticalCoord(const atlas::Field &, const atlas::Field &);
  void setupInterpolation();
  void setupParallelization();
  void setupKernels();
  void setupNormalization();

  // I/O
  void read(const int &);
  void broadcast();
  std::vector<int> writeDef(const int &) const;
  void writeData(const std::vector<int> &) const;

  // Accessors
  double & rh() {return rh_;}
  const double & rh() const {return rh_;}
  const double & rv() const {return rv_;}
  const std::vector<double> & normVertCoord() const {return normVertCoord_;}
  const atlas::FieldSet & norm() const {return norm_;}
  const atlas::FieldSet & normAcc() const {return normAcc_;}

  // Interpolations
  void interpolationTL(const atlas::Field &, atlas::Field &) const;
  void interpolationAD(const atlas::Field &, atlas::Field &) const;

  // Transforms
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

  // Multiply
  void multiplyRedSqrt(const atlas::Field &, atlas::Field &) const;
  void multiplyRedSqrtTrans(const atlas::Field &, atlas::Field &) const;
  void multiplySqrt(const atlas::Field &, atlas::Field &) const;
  void multiplySqrtTrans(const atlas::Field &, atlas::Field &) const;
  void multiply(atlas::Field &) const;

 private:
  // Parameters
  ParametersBase_ params_;

  // Model grid geometry data
  const oops::GeometryData & gdata_;

  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Variable
  std::string myVar_;

  // Model grid
  size_t nx0_;
  size_t ny0_;
  size_t nodes0_;
  size_t nz0_;

  // Reduction factor
  double xRedFac_;
  double yRedFac_;
  double zRedFac_;

  // Reduced grid
  size_t nx_;
  size_t ny_;
  size_t nodes_;
  size_t nz_;
  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;

  // Reduced grid <=> model grid
  size_t rSize_;
  size_t mSize_;
  std::vector<int> rSendCounts_;
  std::vector<int> rSendDispls_;
  std::vector<int> rRecvCounts_;
  std::vector<int> rRecvDispls_;
  std::vector<size_t> sendMapping_;
  std::vector<InterpElement> horInterp_;
  std::vector<InterpElement> verInterp_;

  // Horizontal parallelization
  std::vector<size_t> nxPerTask_;
  std::vector<size_t> nyPerTask_;
  std::vector<size_t> nxStart_;
  std::vector<size_t> nyStart_;
  std::vector<size_t> nxEnd_;
  std::vector<size_t> nyEnd_;

  // Reduced grid <=> rows
  size_t xSize_;
  std::vector<int> xTask_;
  std::vector<int> xOffset_;
  std::vector<int> xSendCounts_;
  std::vector<int> xSendDispls_;
  std::vector<int> xRecvCounts_;
  std::vector<int> xRecvDispls_;
  std::vector<int> xIndex_i_;
  std::vector<int> xIndex_j_;

  // Rows <=> columns
  size_t ySize_;
  std::vector<std::vector<int>> yTask_;
  std::vector<std::vector<int>> yOffset_;
  std::vector<int> ySendCounts_;
  std::vector<int> ySendDispls_;
  std::vector<int> yRecvCounts_;
  std::vector<int> yRecvDispls_;
  std::vector<int> yIndex_i_;
  std::vector<int> yIndex_j_;

  // Convolution
  double rh_;
  double rv_;
  double resol_;
  std::vector<double> normVertCoord_;
  size_t xKernelSize_;
  size_t yKernelSize_;
  size_t zKernelSize_;
  std::vector<double> xKernel_;
  std::vector<double> yKernel_;
  std::vector<double> zKernel_;
  size_t xNormSize_;
  size_t yNormSize_;
  size_t zNormSize_;
  std::vector<double> xNorm_;
  std::vector<double> yNorm_;
  std::vector<double> zNorm_;
  atlas::FieldSet norm_;
  atlas::FieldSet normAcc_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
