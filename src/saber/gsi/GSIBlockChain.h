/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"

#include "saber/gsi/covariance/Covariance.interface.h"

namespace saber {

namespace gsi {

/// GSI covariance block chain with interpolation (optional). Elevated to block
/// chain status because it handles ensemble covariance within.
class SaberGSIBlockChain : public SaberBlockChainBase {
 public:
  template<typename MODEL>
  SaberGSIBlockChain(const oops::Geometry<MODEL> & geom,
                     const oops::Geometry<MODEL> & dualResGeom,
                     const oops::Variables & outerVars,
                     oops::FieldSet4D & fset4dXb,
                     oops::FieldSet4D & fset4dFg,
                     oops::FieldSets & fsetEns,
                     oops::FieldSets & fsetDualResEns,
                     const eckit::LocalConfiguration & covarConf,
                     const eckit::Configuration & conf);
  ~SaberGSIBlockChain();

  /// @brief Randomize the increment according to this B matrix.
  void randomize(oops::FieldSet4D &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(oops::FieldSet4D &) const;
  /// @brief Get this B matrix square-root control vector size.
  size_t ctlVecSize() const {
    throw eckit::NotImplemented(Here());
  }
  /// @brief Multiply the control vector by this B matrix square-root.
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const {
    throw eckit::NotImplemented(Here());
  }
  /// @brief Multiply the increment by this B matrix square-root adjoint.
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const {
    throw eckit::NotImplemented(Here());
  }

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const {return outerVariables_;}

 private:
  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  /// Outer blocks (typically GSI interpolation, but not limited to that)
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;

  // Fortran LinkedList key
  CovarianceKey keySelf_;
  // GSI variables
  oops::Variables centralVars_;
  // GSI grid FunctionSpace
  atlas::FunctionSpace centralFunctionSpace_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberGSIBlockChain::SaberGSIBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResGeom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       oops::FieldSets & fsetEns,
                       oops::FieldSets & fsetDualResEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars) {
  oops::Log::trace() << "SaberGSIBlockChain ctor starting" << std::endl;

  // Check that parallel time decomposition is not used for 4D covariances
  // (currently not supported)
  if (fset4dXb.commTime().size() > 1) {
    throw eckit::NotImplemented("SABER GSI covariance block chain currently does not support "
                                "4DEnVar with parallel decomposition in time");
  }

  // If needed create outer block chain
  if (conf.has("saber outer blocks")) {
    std::vector<SaberOuterBlockParametersWrapper> cmpOuterBlocksParams;
    for (const auto & cmpOuterBlockConf : conf.getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
      cmpOuterBlockParamsWrapper.deserialize(cmpOuterBlockConf);
      cmpOuterBlocksParams.push_back(cmpOuterBlockParamsWrapper);
    }
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVariables_,
                          fset4dXb, fset4dFg, fsetEns, covarConf,
                          cmpOuterBlocksParams);
  }

  // Set outer variables and geometry data for central block (GSI covariance)
  const oops::Variables currentOuterVars = outerBlockChain_ ?
                             outerBlockChain_->innerVars() : outerVariables_;
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geom.generic();

  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Save central function space and variables
  centralFunctionSpace_ = currentOuterGeom.functionSpace();
  centralVars_ = currentOuterVars;

  // Create covariance module
  std::vector<const atlas::field::FieldSetImpl*> fset4dXbptrs(fset4dXb.size());
  std::vector<const atlas::field::FieldSetImpl*> fset4dFgptrs(fset4dFg.size());
  std::vector<const util::DateTime *> timesptrs(fset4dXb.size());
  const std::vector<util::DateTime> & times = fset4dXb.times();
  for (size_t itime = 0; itime < fset4dXb.size(); ++itime) {
    fset4dXbptrs[itime] = fset4dXb[itime].get();
    fset4dFgptrs[itime] = fset4dFg[itime].get();
    timesptrs[itime]    = &times[itime];
  }
  gsi_covariance_create_f90(keySelf_, currentOuterGeom.comm(),
                            saberCentralBlockParams.readParams.value().value(),
                            fset4dXbptrs.size(), fset4dXbptrs.data(), fset4dFgptrs.data(),
                            timesptrs.data());
  // Adjoint test
  // TODO(Anna): this code is similar to the code in CentralBlock::adjoint;
  // the test code need to be generalized so it can be called from different places.
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    // Create random FieldSets
    oops::FieldSet3D fset1_3d = oops::randomFieldSet3D(fset4dXb[0].validTime(),
                                                    geom.generic().comm(),
                                                    geom.generic().functionSpace(),
                                                    outerVars);
    oops::FieldSet3D fset2_3d = oops::randomFieldSet3D(fset4dXb[0].validTime(),
                                                    geom.generic().comm(),
                                                    geom.generic().functionSpace(),
                                                    outerVars);
    oops::FieldSet4D fset1(fset1_3d);
    oops::FieldSet4D fset2(fset2_3d);

    // Copy FieldSets
    oops::FieldSet4D fset1Save(fset1);
    oops::FieldSet4D fset2Save(fset2);

    // Apply forward multiplication only (self-adjointness test)
    this->multiply(fset1);
    this->multiply(fset2);

    // Compute adjoint test
    const double dp1 = fset1.dot_product_with(fset2Save, centralVars_);
    const double dp2 = fset2.dot_product_with(fset1Save, centralVars_);
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                      << localAdjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for block gsi hybrid covariance";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localAdjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Adjoint test failure for block gsi hybrid covariance", Here());
    }
  }

  oops::Log::trace() << "SaberGSIBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace gsi

}  // namespace saber
