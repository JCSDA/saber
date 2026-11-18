/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BiperiodizationImplParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BiperiodizationImplParameters, oops::Parameters)

 public:
  // X-direction inner extension zone
  oops::Parameter<size_t> innerExtNx{"x inner extension", 0, this};

  // Y-direction inner extension zone
  oops::Parameter<size_t> innerExtNy{"y inner extension", 0, this};

  // X-direction outer extension zone
  oops::Parameter<size_t> outerExtNx{"x outer extension", 0, this};

  // Y-direction outer extension zone
  oops::Parameter<size_t> outerExtNy{"y outer extension", 0, this};

  // Inner partitioner (required if outer partitioner is "custom")
  oops::OptionalParameter<std::string> innerPartitioner{"inner partitioner", this};

  // Mixing size
  oops::Parameter<size_t> nmix{"mixing size", 0, this};

  // Mixing scale
  oops::Parameter<double> Lmix{"mixing scale", 1.0, this};

  // Boyd scale
  oops::Parameter<double> Lboyd{"boyd scale", 2.0, this};
};

// -----------------------------------------------------------------------------

class BiperiodizationData {
 public:
  static const std::string classname()
    {return "saber::bifourier::BiperiodizationData";}

  BiperiodizationData(const eckit::mpi::Comm & comm,
                      const atlas::functionspace::StructuredColumns &,
                      const std::vector<int> &,
                      const atlas::functionspace::StructuredColumns &,
                      const size_t &,
                      const size_t &,
                      const size_t &,
                      const size_t &,
                      const size_t &,
                      const double &,
                      const double &,
                      const oops::Variables &);
  ~BiperiodizationData()
    {};

  // Multiply
  void multiply(atlas::FieldSet &) const;

  // Multiply adjoint
  void multiplyAD(atlas::FieldSet &) const;

 private:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Input and output geometries
  atlas::functionspace::StructuredColumns inputFs_;
  std::vector<int> inputPartition_;
  atlas::functionspace::StructuredColumns outputFs_;

  // Variables
  const oops::Variables vars_;

  // Total number of levels (sum of all levels of all active variables)
  size_t nvz_;

  // Counter
  size_t outputJnode_;

  // Biperiodization operations
  size_t localBiperSize_;
  std::vector<size_t> localOutputJnodeVec_;
  std::vector<size_t> localInputJnodeVec_;
  std::vector<double> localWeightVec_;
  size_t commBiperSize_;
  std::vector<size_t> outputJnodeVec_;
  std::vector<size_t> inputJnodeGlbVec_;
  std::vector<size_t> inputTaskVec_;
  std::vector<double> weightVec_;

  // Communication
  size_t sendSize_;
  size_t recvSize_;
  std::vector<int> sendCounts_;
  std::vector<int> sendDispls_;
  std::vector<int> recvCounts_;
  std::vector<int> recvDispls_;

  // Multiply vectors
  std::vector<size_t> sendInputJnode_;
  std::vector<size_t> outputJnodeVecOrdered_;
  std::vector<double> weightVecOrdered_;
  std::vector<size_t> mappingFull2Red_;

  // Private methods

  // Add element to send
  void addBiperElement(const size_t &,
                       const size_t &,
                       const double &);
};

// -----------------------------------------------------------------------------

class BiperiodizationImpl {
 public:
  static const std::string classname()
    {return "saber::bifourier::BiperiodizationImpl";}

  typedef BiperiodizationImplParameters Parameters_;

  BiperiodizationImpl(const oops::GeometryData &,
                      const oops::Variables &,
                      const Parameters_ &);
  ~BiperiodizationImpl()
    {};

  const atlas::FunctionSpace & innerFunctionSpace()
    {return innerFs_;}

  void multiply(atlas::FieldSet & fset) const
    {direct_->multiply(fset);}
  void multiplyAD(atlas::FieldSet & fset) const
    {direct_->multiplyAD(fset);}
  void inverseMultiply(atlas::FieldSet & fset) const
    {inverse_->multiply(fset);}
  void inverseMultiplyAD(atlas::FieldSet & fset) const
    {inverse_->multiplyAD(fset);}

  bool sameFs() const
    {return sameFs_;}

 private:
  // Same grid flag
  bool sameFs_;

  // Inner FunctionSpace
  atlas::functionspace::StructuredColumns innerFs_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Data structure
  std::unique_ptr<BiperiodizationData> direct_;
  std::unique_ptr<BiperiodizationData> inverse_;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
