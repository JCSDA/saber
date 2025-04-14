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

  // Inner partitioner (if outer partitioner is "custom")
  oops::Parameter<std::string> innerPartitioner{"inner partitioner", "checkerboard", this};

  // Mixing size (have an impact if lower than nxExt or nyExt)
  oops::OptionalParameter<size_t> nmix{"mixing size", this};

  // Mixing scale
  oops::Parameter<double> Lmix{"mixing scale", 1.0, this};

  // Boyd scale
  oops::Parameter<double> Lboyd{"boyd scale", 2.0, this};
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

  void multiply(atlas::FieldSet &) const;
  void multiplyAD(atlas::FieldSet &) const;
  void leftInverseMultiply(atlas::FieldSet &) const;

 private:
  // Inner grid
  atlas::StructuredGrid innerGrid_;

  // Inner partition
  std::vector<int> innerPartition_;

  // Inner FunctionSpace
  atlas::functionspace::StructuredColumns innerFs_;

  // Outer geometry data
  const oops::GeometryData & outerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Variables
  const oops::Variables & vars_;

  // Parameters
  Parameters_ params_;

  // Total number of levels (sum of all levels of all active variables)
  size_t nvz_;

  // Biperiodization operations
  size_t localBiperSize_;
  std::vector<size_t> localOuterJnodeVec_;
  std::vector<size_t> localInnerJnodeVec_;
  std::vector<double> localWeightVec_;
  size_t commBiperSize_;
  std::vector<size_t> outerJnodeVec_;
  std::vector<size_t> innerJnodeGlbVec_;
  std::vector<size_t> innerTaskVec_;
  std::vector<double> weightVec_;

  // Communication
  size_t sendSize_;
  size_t recvSize_;
  std::vector<int> sendCounts_;
  std::vector<int> sendDispls_;
  std::vector<int> recvCounts_;
  std::vector<int> recvDispls_;

  // Multiply vectors
  std::vector<size_t> sendInnerJnode_;
  std::vector<size_t> outerJnodeVecOrdered_;
  std::vector<double> weightVecOrdered_;
  std::vector<size_t> mappingFull2Red_;

  // Private methods

  // Add element to send
  void addBiperElement(const size_t &,
                       const size_t &,
                       const size_t &,
                       const double &);
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
