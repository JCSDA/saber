/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"
#include "oops/interface/ModelData.h"
#include "oops/util/Random.h"

#include "saber/oops/SaberBlockChainBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockChain.h"
#include "saber/oops/SaberParametricBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

/// Chain of outer (optional) and an ensemble "block".
class SaberEnsembleBlockChain : public SaberBlockChainBase {
 public:
  template<typename MODEL>
  SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                          const oops::Geometry<MODEL> & dualResGeom,
                          const oops::Variables & outerVars,
                          const atlas::FieldSet & fsetXb,
                          const atlas::FieldSet & fsetFg,
                          const util::DateTime & validTimeOfXbFg,
                          std::vector<atlas::FieldSet> & fsetEns,
                          std::vector<atlas::FieldSet> & dualResolutionFsetEns,
                          const eckit::LocalConfiguration & covarConf,
                          const eckit::Configuration & conf);
  ~SaberEnsembleBlockChain() = default;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(atlas::FieldSet &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(atlas::FieldSet &) const;

 private:
  /// @brief Outer blocks (optional).
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// @brief Localization block chain (optional).
  std::unique_ptr<SaberParametricBlockChain> locBlockChain_;
  /// @brief Ensemble used in the ensemble covariance.
  std::vector<atlas::FieldSet> ensemble_;
  /// @brief Variables used in the ensemble covariance.
  /// TODO(AS): check whether this is needed or can be inferred from ensemble.
  oops::Variables vars_;
  /// @brief Geometry communicator.
  /// TODO(AS): this can be removed once FieldSet4D/FieldSet3D are used.
  const eckit::mpi::Comm & comm_;
  int seed_ = 7;  // For reproducibility
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberEnsembleBlockChain::SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResolutionGeom,
                       const oops::Variables & outerVars,
                       const atlas::FieldSet & fsetXb,
                       const atlas::FieldSet & fsetFg,
                       const util::DateTime & validTimeOfXbFg,
                       // TODO(AS): remove as argument: this should be read inside the
                       // block.
                       std::vector<atlas::FieldSet> & fsetEns,
                       // TODO(AS): remove as argument: this is currently not used (and
                       // when used should be read inside the block.
                       std::vector<atlas::FieldSet> & dualResolutionFsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf)
  : comm_(geom.getComm()) {
  oops::Log::trace() << "SaberEnsembleBlockChain ctor starting" << std::endl;

  // Check that there is an ensemble of at least 2 members.
  if (fsetEns.size() < 2) {
    throw eckit::BadParameter("Ensemble for SaberEnsembleBlockChain has to have at least"
                              " two members.", Here());
  }

  // Create outer blocks if needed
  if (conf.has("saber outer blocks")) {
    std::vector<SaberOuterBlockParametersWrapper> cmpOuterBlocksParams;
    for (const auto & cmpOuterBlockConf : conf.getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
      cmpOuterBlockParamsWrapper.deserialize(cmpOuterBlockConf);
      cmpOuterBlocksParams.push_back(cmpOuterBlockParamsWrapper);
    }
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVars,
                          fsetXb, fsetFg, validTimeOfXbFg, fsetEns, covarConf,
                          cmpOuterBlocksParams);
  }

  // Outer variables and geometry for the ensemble covariance
  const oops::Variables currentOuterVars = outerBlockChain_ ?
                                           outerBlockChain_->innerVars() : outerVars;
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                                     outerBlockChain_->innerGeometryData() : geom.generic();

  // Get parameters:
  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));
  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Get active variables
  const oops::Variables activeVars = getActiveVars(saberCentralBlockParams, currentOuterVars);
  vars_ += activeVars;
  // Check that active variables are present in variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(currentOuterVars.has(var));
  }

  // Ensemble configuration
  eckit::LocalConfiguration ensembleConf
    = covarConf.getSubConfiguration("ensemble configuration");
  // Read inflation field
  eckit::LocalConfiguration centralBlockConf = conf.getSubConfiguration("saber central block");
  const double inflationValue = centralBlockConf.getDouble("inflation value", 1);
  oops::Log::info() << "Info     : Read inflation field" << std::endl;
  atlas::FieldSet inflationField;
  // Read ATLAS inflation file
  if (centralBlockConf.has("inflation field.atlas file")) {
    eckit::LocalConfiguration inflationConf =
                              centralBlockConf.getSubConfiguration("inflation field.atlas file");
    // Read file
    util::readFieldSet(currentOuterGeom.comm(),
                       currentOuterGeom.functionSpace(),
                       activeVars,
                       inflationConf,
                       inflationField);
    // Set name
    inflationField.name() = "inflation";

    // Print FieldSet norm
    oops::Log::test() << "Norm of input parameter inflation: "
                      << util::normFieldSet(inflationField,
                                            activeVars.variables(),
                                            currentOuterGeom.comm())
                      << std::endl;
  }
  // Use model inflation file
  if (centralBlockConf.has("inflation field.model file")) {
    eckit::LocalConfiguration inflationConf =
                              centralBlockConf.getSubConfiguration("inflation field.model file");
    // Copy file
    // Read fieldsets as increments
    // Create increment
    oops::Increment<MODEL> dx(geom, activeVars, validTimeOfXbFg);
    dx.read(inflationConf);
    oops::Log::test() << "Norm of input parameter inflation"
                      << ": " << dx.norm() << std::endl;
    util::copyFieldSet(dx.fieldSet(), inflationField);
  }

  // Apply inflation on ensemble members
  oops::Log::info() << "Info     : Apply inflation on ensemble members" << std::endl;
  for (auto & fset : fsetEns) {
    // Apply local inflation
    if (!inflationField.empty()) {
      util::multiplyFieldSets(fset, inflationField);
    }

    // Apply global inflation
    util::multiplyFieldSet(fset, inflationValue);
  }

  // Ensemble transform
  // For ensemble transform and localization set ensemble size to zero (BUMP needs that)
  // TODO(AS): check if this is used/needed.
  eckit::LocalConfiguration covarConfUpdated(covarConf);
  covarConfUpdated.set("ensemble configuration.ensemble size", 0);
  // Turn off adjoint test for backwards compatibility.
  // TODO(AS): revisit once the way parameters are passed around is refactored.
  covarConfUpdated.set("adjoint test", false);
  const auto & ensTransConf = saberCentralBlockParams.ensembleTransform.value();
  if (ensTransConf != boost::none) {
    oops::Log::info() << "found ens transform: " << *ensTransConf << std::endl;
    // Initialize ensemble transform blockchain
    std::vector<SaberOuterBlockParametersWrapper> ensTransOuterBlocksParams;
    for (const auto & ensTransOuterBlockConf :
                      ensTransConf->getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper ensTransOuterBlockParamsWrapper;
      ensTransOuterBlockParamsWrapper.deserialize(ensTransOuterBlockConf);
      ensTransOuterBlocksParams.push_back(ensTransOuterBlockParamsWrapper);
    }
    std::unique_ptr<SaberOuterBlockChain> ensTransBlockChain =
           std::make_unique<SaberOuterBlockChain>(geom,
             outerVars, fsetXb, fsetFg, validTimeOfXbFg, fsetEns,
             covarConfUpdated, ensTransOuterBlocksParams);

    // Left inverse of ensemble transform on ensemble members
    oops::Log::info() << "Info     : Left inverse of ensemble transform on ensemble members"
                      << std::endl;
    for (auto & fset : fsetEns) {
      ensTransBlockChain->leftInverseMultiply(fset);
    }

    // Add ensemble transform blocks to outer blocks
    // TODO(AS): refactor so there is no need for non-const accessor to outerBlocks
    // in SaberOuterBlockChain.
    oops::Log::info() << "Info     : Add ensemble transform blocks to outer blocks"
                      << std::endl;
    if (outerBlockChain_) {
      std::move(ensTransBlockChain->outerBlocks().begin(),
                ensTransBlockChain->outerBlocks().end(),
                std::back_inserter(outerBlockChain_->outerBlocks()));
      ensTransBlockChain->outerBlocks().erase(ensTransBlockChain->outerBlocks().begin(),
                                              ensTransBlockChain->outerBlocks().end());
    } else {
      outerBlockChain_ = std::move(ensTransBlockChain);
    }
  }

  // Localization
  const auto & locConf = saberCentralBlockParams.localization.value();
  if (locConf != boost::none) {
  // Initialize localization blockchain
     locBlockChain_ = std::make_unique<SaberParametricBlockChain>(geom, dualResolutionGeom,
             outerVars, fsetXb, fsetFg, validTimeOfXbFg, fsetEns, dualResolutionFsetEns,
             covarConfUpdated, *locConf);
  }
  // Direct calibration
  oops::Log::info() << "Info     : Direct calibration" << std::endl;

  // Initialize ensemble
  for (const auto & fset : fsetEns) {
    ensemble_.push_back(fset);
  }

  // Adjoint test
  // TODO(AS): this is now a copy of the test in SaberCentralBlock; needs to be generalized.
  // (Perhaps the adjoint[s] test can be moved to SaberBlockChainBase.
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));
    // Create random FieldSets
    atlas::FieldSet fset1 =  util::createRandomFieldSet(currentOuterGeom.comm(),
                                                        currentOuterGeom.functionSpace(),
                                                        activeVars);
    atlas::FieldSet fset2 =  util::createRandomFieldSet(currentOuterGeom.comm(),
                                                        currentOuterGeom.functionSpace(),
                                                        activeVars);
    // Copy FieldSets
    atlas::FieldSet fset1Save = util::copyFieldSet(fset1);
    atlas::FieldSet fset2Save = util::copyFieldSet(fset2);

    // Apply forward multiplication only (self-adjointness test)
    // TODO(AS): need to change this to only call it for the ensemble part, not outer blocks!
    this->multiply(fset1);
    this->multiply(fset2);

    // Compute adjoint test
    const double dp1 = util::dotProductFieldSets(fset1, fset2Save,
                                                 activeVars.variables(), currentOuterGeom.comm());
    const double dp2 = util::dotProductFieldSets(fset2, fset1Save,
                                                 activeVars.variables(), currentOuterGeom.comm());
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                      << localAdjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for block Ensemble";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localAdjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      ABORT("Adjoint test failure for block Ensemble");
    }
  }
  oops::Log::trace() << "SaberEnsembleBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
