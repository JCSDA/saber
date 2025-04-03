/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"

#include "saber/fastlam/FastLAMParametersBase.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class LayerBase : public util::Printable,
                  private boost::noncopyable {
 public:
  static const std::string classname() {return "saber::fastlam::Layer";}

  // Constructor
  LayerBase(const FastLAMParametersBase & params,
            const eckit::LocalConfiguration & fieldsMetaData,
            const oops::GeometryData & gdata,
            const std::string & myGroup,
            const std::vector<std::string> & myVars,
            const size_t & nx0,
            const size_t & ny0,
            const size_t & nz0) :
    params_(params),
    fieldsMetaData_(fieldsMetaData),
    gdata_(gdata),
    comm_(gdata_.comm()),
    myrank_(comm_.rank()),
    myGroup_(myGroup),
    myVars_(myVars),
    nx0_(nx0),
    ny0_(ny0),
    mSize_(gdata_.functionSpace().ghost().shape(0)),
    nz0_(nz0),
    parallelization_(params.parallelization.value()) {}
  virtual ~LayerBase() {}

  // Setup
  virtual void setupParallelization() = 0;
  virtual void extractConvolution(const size_t &,
                                  const size_t &,
                                  std::vector<double> &,
                                  std::vector<double> &) = 0;

  // Multiply square-root and adjoint
  virtual size_t ctlVecSize() const = 0;
  virtual void multiplySqrt(const atlas::Field &,
                            atlas::Field &,
                            const size_t &) const = 0;
  virtual void multiplySqrtTrans(const atlas::Field &,
                                 atlas::Field &,
                                 const size_t &) const = 0;

  // Non-virtual methods

  // Setups
  void setupVerticalCoord(const atlas::Field &,
                          const atlas::Field &);
  void setupInterpolation();
  void testInterpolation(const std::vector<double> &) const;
  void setupKernels();
  void setupNormalization();

  // I/O
  void read(const int &);
  void broadcast();
  std::array<int, 8> writeDef(const int &) const;
  void writeData(const std::array<int, 8> &) const;

  // Accessors
  double & rh() {return rh_;}
  const double & rh() const {return rh_;}
  double & rv() {return rv_;}
  const double & rv() const {return rv_;}
  double & resol() {return resol_;}
  const double & resol() const {return resol_;}
  double & rfh() {return rfh_;}
  const double & rfh() const {return rfh_;}
  double & rfv() {return rfv_;}
  const double & rfv() const {return rfv_;}
  const std::vector<double> & normVertCoord() const {return normVertCoord_;}
  const atlas::FieldSet & norm() const {return norm_;}
  const atlas::FieldSet & normAcc() const {return normAcc_;}

 protected:
  // Interpolations
  void interpolationTL(const atlas::Field &, atlas::Field &) const;
  void interpolationAD(const atlas::Field &, atlas::Field &) const;

  // Parameters
  FastLAMParametersBase params_;
  const eckit::LocalConfiguration fieldsMetaData_;

  // Model grid geometry data
  const oops::GeometryData & gdata_;

  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Group and variables
  std::string myGroup_;
  std::vector<std::string> myVars_;

  // Model grid
  size_t nx0_;
  size_t ny0_;
  size_t mSize_;
  size_t nz0_;

  // Resolution and reduction factors
  double resol_;
  double rfh_;
  double rfv_;

  // Convolution
  double rh_;
  double rv_;
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

  // Reduction factor
  double xRedFac_;
  double yRedFac_;
  double zRedFac_;

  // Reduced grid
  size_t nx_;
  size_t ny_;
  size_t rSize_;
  size_t nz_;
  std::vector<int> mpiTask_;
  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;

  // Reduced grid <=> model grid (interpolation)
  bool noInterp_;
  size_t rSendSize_;
  size_t mRecvSize_;
  std::vector<int> rSendCounts_;
  std::vector<int> rSendDispls_;
  std::vector<int> mRecvCounts_;
  std::vector<int> mRecvDispls_;
  std::vector<size_t> rSendMapping_;
  std::vector<std::string> horType_;
  std::vector<std::array<size_t, 4>> horStencil_;
  std::vector<std::array<double, 4>> horWeights_;
  std::vector<size_t> horStencilSize_;
  std::vector<size_t> horIndexX_;
  std::vector<size_t> horIndexY_;
  std::vector<std::array<size_t, 2>> verStencil_;
  std::vector<std::array<double, 2>> verWeights_;
  std::vector<size_t> verStencilSize_;
  std::vector<size_t> verIndex_;

 private:
  // Parallelization mode
  std::string parallelization_;
  virtual void print(std::ostream &) const = 0;
  void binarySearch(const std::vector<int> &,
                    const std::vector<size_t> &,
                    const size_t &,
                    int &);
};

// -----------------------------------------------------------------------------

class LayerFactory;

// -----------------------------------------------------------------------------

class LayerFactory {
 public:
  static std::unique_ptr<LayerBase> create(const FastLAMParametersBase &,
                                           const eckit::LocalConfiguration &,
                                           const oops::GeometryData &,
                                           const std::string &,
                                           const std::vector<std::string> &,
                                           const size_t &,
                                           const size_t &,
                                           const size_t &);

  virtual ~LayerFactory() = default;

 protected:
  explicit LayerFactory(const std::string &name);

 private:
  virtual std::unique_ptr<LayerBase> make(const FastLAMParametersBase &,
                                          const eckit::LocalConfiguration &,
                                          const oops::GeometryData &,
                                          const std::string &,
                                          const std::vector<std::string> &,
                                          const size_t &,
                                          const size_t &,
                                          const size_t &) = 0;

  static std::map < std::string, LayerFactory * > & getMakers() {
    static std::map < std::string, LayerFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LayerMaker : public LayerFactory {
  std::unique_ptr<LayerBase> make(const FastLAMParametersBase & params,
                                  const eckit::LocalConfiguration & fieldsMetaData,
                                  const oops::GeometryData & gdata,
                                  const std::string & myGroup,
                                  const std::vector<std::string> & myVars,
                                  const size_t & nx0,
                                  const size_t & ny0,
                                  const size_t & nz0) override {
    return std::make_unique<T>(params, fieldsMetaData, gdata, myGroup, myVars, nx0, ny0, nz0);
  }

 public:
  explicit LayerMaker(const std::string & name) : LayerFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
