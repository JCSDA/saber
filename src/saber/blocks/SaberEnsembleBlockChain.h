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

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/interface/ModelData.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

/// Chain of outer (optional) and an ensemble "block".
class SaberEnsembleBlockChain : public SaberBlockChainBase {
 public:
  template<typename MODEL>
  SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                          const oops::Geometry<MODEL> & dualResGeom,
                          const oops::Variables & outerVars,
                          oops::FieldSet4D & fset4dXb,
                          oops::FieldSet4D & fset4dFg,
                          oops::FieldSets & fsetEns,
                          oops::FieldSets & fsetDualResEns,
                          const eckit::LocalConfiguration & covarConf,
                          const eckit::Configuration & conf);
  ~SaberEnsembleBlockChain() = default;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(oops::FieldSet4D &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(oops::FieldSet4D &) const;
  /// @brief Get this B matrix square-root control vector size.
  size_t ctlVecSize() const {return ctlVecSize_;}
  /// @brief Multiply the control vector by this B matrix square-root.
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const;
  /// @brief Multiply the increment by this B matrix square-root adjoint.
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const;

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const {return outerVariables_;}

 private:
  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  /// @brief Outer blocks (optional).
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// @brief Localization block chain (optional).
  std::unique_ptr<SaberParametricBlockChain> locBlockChain_;
  /// @brief Ensemble used in the ensemble covariance.
  oops::FieldSets ensemble_;
  /// @brief Control vector size.
  size_t ctlVecSize_;
  /// @brief Variables used in the ensemble covariance.
  /// TODO(AS): check whether this is needed or can be inferred from ensemble.
  oops::Variables vars_;
  int seed_ = 7;  // For reproducibility
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberEnsembleBlockChain::SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResGeom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       // TODO(AS): remove as argument: this should be read inside the
                       // block.
                       oops::FieldSets & fsetEns,
                       // TODO(AS): remove as argument: this is currently not used (and
                       // when used should be read inside the block.
                       oops::FieldSets & fsetDualResEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars),
    ensemble_(fsetEns), ctlVecSize_(0) {
  oops::Log::trace() << "SaberEnsembleBlockChain ctor starting" << std::endl;

  // Check that there is an ensemble of at least 2 members.
  if (ensemble_.ens_size() < 2) {
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
                          fset4dXb, fset4dFg, ensemble_, covarConf,
                          cmpOuterBlocksParams);
  }

  // Outer variables and geometry for the ensemble covariance
  oops::Variables currentOuterVars = outerBlockChain_ ?
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
  for (const auto & var : activeVars) {
    if (!currentOuterVars.has(var)) {
      throw eckit::UserError("Active variable " + var.name() + " is not present in "
                             "outer variables", Here());
    }
  }

  // Ensemble configuration
  eckit::LocalConfiguration ensembleConf
    = covarConf.getSubConfiguration("ensemble configuration");

  // Check consistency if ensemble was read on non-MODEL geometry
  if (ensembleConf.has("ensemble pert on other geometry")) {
    const auto & currentFspace = currentOuterGeom.functionSpace();
    const auto & ensFspace = ensemble_[0].fieldSet()[0].functionspace();
    ASSERT(ensFspace.type() == currentFspace.type());
    ASSERT(util::getGridUid(ensFspace).compare(
               util::getGridUid(currentFspace)) == 0);
  }

  // Read inflation field
  eckit::LocalConfiguration centralBlockConf = conf.getSubConfiguration("saber central block");
  const double inflationValue = centralBlockConf.getDouble("inflation value", 1);
  oops::Log::info() << "Info     : Read inflation field" << std::endl;
  oops::FieldSet3D inflationField(fset4dXb[0].validTime(), geom.getComm());
  // Read ATLAS inflation file
  if (centralBlockConf.has("inflation field.atlas file")) {
    eckit::LocalConfiguration inflationConf =
                              centralBlockConf.getSubConfiguration("inflation field.atlas file");
    // Read file
    inflationField.read(currentOuterGeom.functionSpace(),
                        activeVars,
                        inflationConf);
    // Set name
    inflationField.name() = "inflation";

    // Print FieldSet norm
    oops::Log::test() << "Norm of input parameter inflation: "
                      << inflationField.norm(activeVars) << std::endl;
  }
  // Use model inflation file
  if (centralBlockConf.has("inflation field.model file")) {
    eckit::LocalConfiguration inflationConf =
                              centralBlockConf.getSubConfiguration("inflation field.model file");
    // Copy file
    // Read fieldsets as increments
    // Create increment
    oops::Increment<MODEL> dx(geom, activeVars, fset4dXb[0].validTime());
    dx.read(inflationConf);
    oops::Log::test() << "Norm of input parameter inflation"
                      << ": " << dx.norm() << std::endl;
    inflationField.deepCopy(dx.fieldSet());
  }

  // Apply inflation on ensemble members
  oops::Log::info() << "Info     : Apply inflation on ensemble members" << std::endl;
  // Apply local inflation
  if (!inflationField.empty()) {
    ensemble_ *= inflationField;
  }
  // Apply global inflation
  ensemble_ *= inflationValue;

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
             currentOuterVars, fset4dXb, fset4dFg, ensemble_,
             covarConfUpdated, ensTransOuterBlocksParams);

    // Left inverse of ensemble transform on ensemble members
    oops::Log::info() << "Info     : Left inverse of ensemble transform on ensemble members"
                      << std::endl;
    for (size_t itime = 0; itime < ensemble_.local_time_size(); ++itime) {
      for (size_t iens = 0; iens < ensemble_.local_ens_size(); ++iens) {
        ensTransBlockChain->leftInverseMultiply(ensemble_(itime, iens));
      }
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

    // Update outer variables
    currentOuterVars = outerBlockChain_->innerVars();
  }

  const oops::GeometryData & localizationOuterGeomData = outerBlockChain_ ?
              outerBlockChain_->innerGeometryData() : geom.generic();

  // Localization
  const auto & locConf = saberCentralBlockParams.localization.value();
  if (locConf != boost::none) {
    // The localization is a parametric block chain constructed with the same geometry
    // as the ensemble block chain by default. If the outer blocks or transform block
    // chain include a change of geometries, we need to use build the localization from
    // the current geometryData, using a generic constructor of the parametric block chain.

    // Check consistency of `geom` and the current geometry
    const auto & currentFspace = localizationOuterGeomData.functionSpace();
    if (util::getGridUid(geom.functionSpace()) != util::getGridUid(currentFspace)) {
      oops::Log::info() << "Info     : Localization and ensemble are on different "
                           "functionSpaces, building localization with generic "
                           "constructor" << std::endl;
      // Note QUENCH could just build another geometry here and use the standard
      // constructor, but other models usually don't have this ability to create a
      // Geometry on any mesh.
      locBlockChain_ = std::make_unique<SaberParametricBlockChain>(localizationOuterGeomData,
                                                                   currentOuterVars,
                                                                   fset4dXb,
                                                                   fset4dFg,
                                                                   covarConfUpdated,
                                                                   *locConf);
    } else {
      oops::Log::info() << "Info     : Localization and ensemble are on same "
                           "functionSpaces, building localization with standard "
                           "constructor" << std::endl;
      locBlockChain_ = std::make_unique<SaberParametricBlockChain>(geom,
                                                                   dualResGeom,
                                                                   currentOuterVars,
                                                                   fset4dXb,
                                                                   fset4dFg,
                                                                   ensemble_,
                                                                   fsetDualResEns,
                                                                   covarConfUpdated,
                                                                   *locConf);
    }
  }
  // Direct calibration
  oops::Log::info() << "Info     : Direct calibration" << std::endl;

  // Get control vector size
  if (locBlockChain_) {
    // With localization
    ctlVecSize_ = ensemble_.ens_size()*locBlockChain_->ctlVecSize();
  } else {
    // Without localization
    ctlVecSize_ = ensemble_.ens_size();
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
    oops::FieldSet4D fset4d1(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d1.size(); ++jtime) {
      fset4d1[jtime].randomInit(currentOuterGeom.functionSpace(), activeVars);
    }
    oops::FieldSet4D fset4d2(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d2.size(); ++jtime) {
      fset4d2[jtime].randomInit(currentOuterGeom.functionSpace(), activeVars);
    }
    // Copy FieldSets
    const oops::FieldSet4D fset4d1Save = oops::copyFieldSet4D(fset4d1);
    const oops::FieldSet4D fset4d2Save = oops::copyFieldSet4D(fset4d2);

    // Apply forward multiplication only (self-adjointness test)
    // TODO(AS): need to change this to only call it for the ensemble part, not outer blocks!
    this->multiply(fset4d1);
    this->multiply(fset4d2);

    // Compute adjoint test
    const double dp1 = fset4d1.dot_product_with(fset4d2Save, activeVars);
    const double dp2 = fset4d2.dot_product_with(fset4d1Save, activeVars);
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                      << localAdjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for block Ensemble";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localAdjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Adjoint test failure for block Ensemble", Here());
    }
  }

  // Square-root test
  if (covarConf.getBool("square-root test")) {
    // Get tolerance
    const double localSqrtTolerance =
      saberCentralBlockParams.sqrtTolerance.value().get_value_or(
      covarConf.getDouble("square-root tolerance"));

    // Create FieldSet
    oops::FieldSet4D fset4d(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      fset4d[jtime].randomInit(outerFunctionSpace_, outerVariables_);
    }

    // Copy FieldSet
    oops::FieldSet4D fset4dSave = oops::copyFieldSet4D(fset4d);

    // Create control vector
    oops::Log::info() << "Control vector size for block Ensemble: "
                      << ctlVecSize() << std::endl;
    atlas::Field ctlVec = atlas::Field("genericCtlVec",
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::make_shape(this->ctlVecSize()));
    size_t seed = 7;  // To avoid impact on future random generator calls
    util::NormalDistribution<double> dist(this->ctlVecSize(), 0.0, 1.0, seed);
    std::vector<double> randVec;
    for (size_t jnode = 0; jnode < this->ctlVecSize(); ++jnode) {
      randVec.push_back(dist[jnode]);
    }
    if (!locBlockChain_) {
      currentOuterGeom.comm().broadcast(randVec, 0);
    }
    auto view = atlas::array::make_view<double, 1>(ctlVec);
    for (size_t jnode = 0; jnode < this->ctlVecSize(); ++jnode) {
      view(jnode) = randVec[jnode];
    }

    // Copy control vector
    atlas::Field ctlVecSave = atlas::Field("genericCtlVec",
                                           atlas::array::make_datatype<double>(),
                                           atlas::array::make_shape(this->ctlVecSize()));
    auto viewSave = atlas::array::make_view<double, 1>(ctlVecSave);
    viewSave.assign(view);

    // Apply square-root multiplication
    this->multiplySqrt(ctlVecSave, fset4d, 0);

    // Apply square-root adjoint multiplication
    this->multiplySqrtAD(fset4dSave, ctlVec, 0);

    // Compute adjoint test
    const double dp1 = fset4d.dot_product_with(fset4dSave, activeVars);
    double dp2 = 0.0;
    for (size_t jnode = 0; jnode < this->ctlVecSize(); ++jnode) {
      dp2 += view(jnode)*viewSave(jnode);
    }
    if (locBlockChain_) {
      currentOuterGeom.comm().allReduceInPlace(dp2, eckit::mpi::sum());
      fset4d.commTime().allReduceInPlace(dp2, eckit::mpi::sum());
    }
    oops::Log::info() << std::setprecision(16) << "Info     : Square-root test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : square-root tolerance = "
                      << localSqrtTolerance << std::endl;
    const bool adjComparison = (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localSqrtTolerance);

    // Apply square-root multiplication
    this->multiplySqrt(ctlVec, fset4d, 0);

    // Apply full multiplication
    this->multiply(fset4dSave);

    // Check that the fieldsets are similar within tolerance
    bool sqrtComparison = true;
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      sqrtComparison = sqrtComparison && fset4d[jtime].compare_with(fset4dSave[jtime],
        localSqrtTolerance, util::ToleranceType::relative);
    }
    if (sqrtComparison) {
      oops::Log::info() << "Info     : Square-root test passed: U U^t x == B x" << std::endl;
    } else {
      oops::Log::info() << "Info     : Square-root test failed: U U^t x != B x" << std::endl;
    }

    // Print results
    oops::Log::test() << "Square-root test for block Ensemble";
    if (adjComparison && sqrtComparison) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Square-root test failure for block Ensemble", Here());
    }
  }

  oops::Log::trace() << "SaberEnsembleBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
