/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberCentralBlockBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"
#include "oops/util/Random.h"

#include "saber/blocks/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberCentralBlockFactory::SaberCentralBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberCentralBlockFactory."
                       << std::endl;
    throw eckit::Exception("Element already registered in saber::SaberCentralBlockFactory.",
      Here());
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberCentralBlockBase> SaberCentralBlockFactory::create(
  const oops::GeometryData & geometryData,
  const oops::Variables & vars,
  const eckit::Configuration & covarConf,
  const SaberBlockParametersBase & params,
  const oops::FieldSet3D & xb,
  const oops::FieldSet3D & fg) {
  oops::Log::trace() << "SaberCentralBlockBase::create starting" << std::endl;
  const std::string id = params.saberBlockName;
  typename std::map<std::string, SaberCentralBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberCentralBlockFactory." << std::endl;
    throw eckit::UserError("Element does not exist in saber::SaberCentralBlockFactory.", Here());
  }
  std::unique_ptr<SaberCentralBlockBase> ptr =
    jsb->second->make(geometryData, vars, covarConf,
                      params, xb, fg);
  oops::Log::trace() << "SaberCentralBlockBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberBlockParametersBase>
SaberCentralBlockFactory::createParameters(const std::string &name) {
  typename std::map<std::string, SaberCentralBlockFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in saber::SaberCentralBlockFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

void SaberCentralBlockBase::adjointTest(const oops::GeometryData & geometryData,
                                        const oops::Variables & vars,
                                        const double & adjointTolerance) const {
  oops::Log::trace() << "SaberCentralBlockBase::adjointTest starting" << std::endl;

  // Create random FieldSets
  oops::FieldSet3D fset1 = oops::randomFieldSet3D(validTime_,
                                                  geometryData.comm(),
                                                  geometryData.functionSpace(),
                                                  vars);
  oops::FieldSet3D fset2 = oops::randomFieldSet3D(validTime_,
                                                  geometryData.comm(),
                                                  geometryData.functionSpace(),
                                                  vars);

  // Copy FieldSets
  oops::FieldSet3D fset1Save(fset1);
  oops::FieldSet3D fset2Save(fset2);

  // Apply forward multiplication only (self-adjointness test)
  this->multiply(fset1);
  this->multiply(fset2);

  // Compute adjoint test
  const double dp1 = fset1.dot_product_with(fset2Save, vars);
  const double dp2 = fset2.dot_product_with(fset1Save, vars);
  oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: (Ax)^t y = " << dp1
                    << ": x^t (Ay) = " << dp2 << " : adjoint tolerance = "
                    << adjointTolerance << std::endl;
  oops::Log::test() << "Adjoint test for block " << this->blockName();
  if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < adjointTolerance) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    throw eckit::Exception("Adjoint test failure for block " + this->blockName(), Here());
  }

  oops::Log::trace() << "SaberCentralBlockBase::adjointTest done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlockBase::sqrtTest(const oops::GeometryData & geometryData,
                                     const oops::Variables & vars,
                                     const double & sqrtTolerance) const {
  oops::Log::trace() << "SaberOuterBlockBase::sqrtTest starting" << std::endl;

  // Square-root test

  // Create FieldSet
  oops::FieldSet3D fset = oops::randomFieldSet3D(validTime_,
                                                 geometryData.comm(),
                                                 geometryData.functionSpace(),
                                                 vars);

  // Copy FieldSet
  oops::FieldSet3D fsetSave(fset);

  // Create control vector
  oops::Log::info() << "Control vector size for block " << this->blockName() << ": "
                      << ctlVecSize() << std::endl;
  atlas::Field ctlVec = atlas::Field("genericCtlVec",
                                     atlas::array::make_datatype<double>(),
                                     atlas::array::make_shape(ctlVecSize()));
  size_t seed = 7;  // To avoid impact on future random generator calls
  util::NormalDistribution<double> dist(ctlVecSize(), 0.0, 1.0, seed);
  auto view = atlas::array::make_view<double, 1>(ctlVec);
  for (size_t jnode = 0; jnode < ctlVecSize(); ++jnode) {
    view(jnode) = dist[jnode];
  }

  // Copy control vector
  atlas::Field ctlVecSave = atlas::Field("genericCtlVec",
                                         atlas::array::make_datatype<double>(),
                                         atlas::array::make_shape(ctlVecSize()));
  auto viewSave = atlas::array::make_view<double, 1>(ctlVecSave);
  viewSave.assign(view);

  // Apply square-root multiplication
  this->multiplySqrt(ctlVecSave, fset, 0);

  // Apply square-root adjoint multiplication
  this->multiplySqrtAD(fsetSave, ctlVec, 0);

  // Compute adjoint test
  const double dp1 = fset.dot_product_with(fsetSave, vars);
  double dp2 = 0.0;
  for (size_t jnode = 0; jnode < ctlVecSize(); ++jnode) {
    dp2 += view(jnode)*viewSave(jnode);
  }
  geometryData.comm().allReduceInPlace(dp2, eckit::mpi::sum());
  oops::Log::info() << std::setprecision(16) << "Info     : Square-root test: y^t (Ux) = " << dp1
                    << ": x^t (U^t y) = " << dp2 << " : square-root tolerance = "
                    << sqrtTolerance << std::endl;
  const bool adjComparison = (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < sqrtTolerance);

  // Apply square-root multiplication
  this->multiplySqrt(ctlVec, fset, 0);

  // Apply full multiplication
  this->multiply(fsetSave);

  // Check that the fieldsets are similar within tolerance
  const bool sqrtComparison = fset.compare_with(fsetSave, sqrtTolerance,
                                                util::ToleranceType::relative);
  if (sqrtComparison) {
    oops::Log::info() << "Info     : Square-root test passed: U U^t x == B x" << std::endl;
  } else {
    oops::Log::info() << "Info     : Square-root test failed: U U^t x != B x" << std::endl;
  }

  // Print results
  oops::Log::test() << "Square-root test for block " << this->blockName();
  if (adjComparison && sqrtComparison) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    throw eckit::Exception("Square-root test failure for block " + this->blockName(), Here());
  }

  oops::Log::trace() << "SaberOuterBlockBase::sqrtTest done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
