/*
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberEnsembleBlockChain.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef oops::Increment4D<MODEL>                             Increment4D_;
  typedef oops::State4D<MODEL>                                 State4D_;
  typedef typename oops::Increment<MODEL>::WriteParameters_    WriteParameters_;

 public:
  typedef ErrorCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &, const Parameters_ &,
                  const State4D_ &, const State4D_ &);
  virtual ~ErrorCovariance();

  void multiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {this->doMultiply(dxi, dxo);}

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;

  /// Chain of outer blocks applied to all components of hybrid covariances.
  /// Not initialized for non-hybrid covariances.
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// Vector of hybrid B components (one element for non-hybrid case).
  std::vector<std::unique_ptr<SaberBlockChainBase>> hybridBlockChain_;
  /// Vector of scalar weights for hybrid B components (one element, equal to
  /// 1.0 for non-hybrid case).
  std::vector<double> hybridScalarWeightSqrt_;
  /// Vector of field weights for hybrid B components (one element, empty
  /// fieldset for non-hybrid case).
  std::vector<atlas::FieldSet> hybridFieldWeightSqrt_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const Parameters_ & params,
                                        const State4D_ & xb,
                                        const State4D_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, params, xb, fg)
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Local copy of background and first guess that can undergo interpolation
  atlas::FieldSet fsetXb = util::copyFieldSet(xb[0].fieldSet());
  atlas::FieldSet fsetFg = util::copyFieldSet(fg[0].fieldSet());
  ASSERT(xb[0].validTime() == fg[0].validTime());
  const util::DateTime validTimeOfXbFg = xb[0].validTime();

  // Extend backgroud and first guess with extra fields
  // TODO(Benjamin, Marek, Mayeul, ?)

  // Initialize outer variables
  const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
  oops::Variables outerVars(incVars.variables());
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    outerVars.addMetaData(outerVars[i], "levels", vlevs[i]);
  }

  // Create covariance configuration
  eckit::LocalConfiguration covarConf;
  covarConf.set("adjoint test", params.adjointTest.value());
  covarConf.set("adjoint tolerance", params.adjointTolerance.value());
  covarConf.set("inverse test", params.inverseTest.value());
  covarConf.set("inverse tolerance", params.inverseTolerance.value());
  covarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = params.iterativeEnsembleLoading.value();

  // Initialize ensembles as vector of FieldSets
  std::vector<atlas::FieldSet> fsetEns;
  // Read ensemble (for non-iterative ensemble loading)
  eckit::LocalConfiguration ensembleConf = readEnsemble(geom,
                                                        outerVars,
                                                        xb[0],
                                                        fg[0],
                                                        params.toConfiguration(),
                                                        iterativeEnsembleLoading,
                                                        fsetEns);
  covarConf.set("ensemble configuration", ensembleConf);
  // Read dual resolution ensemble if needed
  const auto & dualResolutionParams = params.dualResolutionParams.value();
  const oops::Geometry<MODEL> * dualResolutionGeom = &geom;
  std::vector<atlas::FieldSet> dualResolutionFsetEns;
  if (dualResolutionParams != boost::none) {
    const auto & dualResolutionGeomConf = dualResolutionParams->geometry.value();
    if (dualResolutionGeomConf != boost::none) {
      // Create dualResolution geometry
      typename oops::Geometry<MODEL>::Parameters_ dualResolutionGeometryParams;
      dualResolutionGeometryParams.deserialize(*dualResolutionGeomConf);
      dualResolutionGeom = new oops::Geometry<MODEL>(dualResolutionGeometryParams, geom.getComm());
    }
    // Background and first guess at dual resolution geometry
    oops::State<MODEL> xbDualResolution(*dualResolutionGeom, xb[0]);
    oops::State<MODEL> fgDualResolution(*dualResolutionGeom, fg[0]);
    // Read dual resolution ensemble
    eckit::LocalConfiguration dualResolutionEnsembleConf
      = readEnsemble(*dualResolutionGeom,
                     outerVars,
                     xbDualResolution,
                     fgDualResolution,
                     dualResolutionParams->toConfiguration(),
                     iterativeEnsembleLoading,
                     dualResolutionFsetEns);

    // Add dual resolution ensemble configuration
    covarConf.set("dual resolution ensemble configuration", dualResolutionEnsembleConf);
  }

  // Add ensemble output
  const auto & outputEnsemble = params.outputEnsemble.value();
  if (outputEnsemble != boost::none) {
    covarConf.set("output ensemble", outputEnsemble->toConfiguration());
  }

  const SaberBlockParametersBase & saberCentralBlockParams =
    params.saberCentralBlockParams.value().saberCentralBlockParameters;
  // Build covariance blocks: hybrid covariance case
  if (saberCentralBlockParams.saberBlockName.value() == "Hybrid") {
    // Build common (for all hybrid components) outer blocks if they exist
    const auto & saberOuterBlocksParams = params.saberOuterBlocksParams.value();
    if (saberOuterBlocksParams != boost::none) {
      outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom,
                       outerVars,
                       fsetXb,
                       fsetFg,
                       validTimeOfXbFg,
                       fsetEns,
                       covarConf,
                       *saberOuterBlocksParams);
      outerVars = outerBlockChain_->innerVars();
    }

    // Hybrid central block
    eckit::LocalConfiguration hybridConf = saberCentralBlockParams.toConfiguration();

    // Create block geometry (needed for ensemble reading)
    const oops::Geometry<MODEL> * hybridGeom = &geom;
    if (hybridConf.has("geometry")) {
      hybridGeom = new oops::Geometry<MODEL>(hybridConf.getSubConfiguration("geometry"),
        geom.getComm());
    }

    // Loop over components
    for (const auto & cmp : hybridConf.getSubConfigurations("components")) {
      // Initialize component outer variables
      // TODO(AS): this should be either outerVars or outerBlockChain_->innerVars();
      oops::Variables cmpOuterVars(outerVars);

      // Set weight
      eckit::LocalConfiguration weightConf = cmp.getSubConfiguration("weight");
      // Scalar weight
      hybridScalarWeightSqrt_.push_back(std::sqrt(weightConf.getDouble("value", 1.0)));
      // File-base weight
      atlas::FieldSet weightFset;
      if (weightConf.has("file")) {
        // File-base weight
        readHybridWeight(*hybridGeom,
                         outerVars,
                         validTimeOfXbFg,
                         weightConf.getSubConfiguration("file"),
                         weightFset);
        util::sqrtFieldSet(weightFset);
      }
      hybridFieldWeightSqrt_.push_back(weightFset);

      // Set covariance
      eckit::LocalConfiguration cmpConf = cmp.getSubConfiguration("covariance");

      // Initialize ensembles as vector of FieldSets
      std::vector<atlas::FieldSet> cmpFsetEns;
      // Read ensemble
      eckit::LocalConfiguration cmpEnsembleConf
         = readEnsemble(*hybridGeom,
                        cmpOuterVars,
                        xb[0],
                        fg[0],
                        cmpConf,
                        params.iterativeEnsembleLoading.value(),
                        cmpFsetEns);

      // Create internal configuration
      eckit::LocalConfiguration cmpCovarConf;
      cmpCovarConf.set("ensemble configuration", cmpEnsembleConf);
      cmpCovarConf.set("adjoint test", params.adjointTest.value());
      cmpCovarConf.set("adjoint tolerance", params.adjointTolerance.value());
      cmpCovarConf.set("inverse test", params.inverseTest.value());
      cmpCovarConf.set("inverse tolerance", params.inverseTolerance.value());
      cmpCovarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());

      SaberCentralBlockParametersWrapper cmpCentralBlockParamsWrapper;
      cmpCentralBlockParamsWrapper.deserialize(cmpConf.getSubConfiguration("saber central block"));
      const auto & centralBlockParams =
                   cmpCentralBlockParamsWrapper.saberCentralBlockParameters.value();
      // TODO(AS): move construction of BlockChain to factory method or function
      if (centralBlockParams.saberBlockName.value() == "Ensemble") {
        hybridBlockChain_.push_back(std::make_unique<SaberEnsembleBlockChain>(*hybridGeom,
                          *dualResolutionGeom,
                          cmpOuterVars,
                          fsetXb,
                          fsetFg,
                          validTimeOfXbFg,
                          cmpFsetEns,
                          dualResolutionFsetEns,
                          cmpCovarConf,
                          cmpConf));
      } else {
        hybridBlockChain_.push_back(std::make_unique<SaberParametricBlockChain>(*hybridGeom,
                          *dualResolutionGeom,
                          cmpOuterVars,
                          fsetXb,
                          fsetFg,
                          validTimeOfXbFg,
                          cmpFsetEns,
                          dualResolutionFsetEns,
                          cmpCovarConf,
                          cmpConf));
      }
    }
  } else {
    // Non-hybrid covariance: single block chain
    if (saberCentralBlockParams.saberBlockName.value() == "Ensemble") {
      hybridBlockChain_.push_back(std::make_unique<SaberEnsembleBlockChain>(geom,
                          *dualResolutionGeom,
                          outerVars,
                          fsetXb,
                          fsetFg,
                          validTimeOfXbFg,
                          fsetEns,
                          dualResolutionFsetEns,
                          covarConf,
                          params.toConfiguration()));
    } else {
      hybridBlockChain_.push_back(std::make_unique<SaberParametricBlockChain>(geom,
                          *dualResolutionGeom,
                          outerVars,
                          fsetXb,
                          fsetFg,
                          validTimeOfXbFg,
                          fsetEns,
                          dualResolutionFsetEns,
                          covarConf,
                          params.toConfiguration()));
    }

    // Set weights
    hybridScalarWeightSqrt_.push_back(1.0);
    // File-base weight
    atlas::FieldSet weightFset;
    hybridFieldWeightSqrt_.push_back(weightFset);
  }

  oops::Log::trace() << "ErrorCovariance::ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance");
  oops::Log::trace() << "ErrorCovariance<MODEL>::~ErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doRandomize(Increment4D_ & dx) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  // SABER block chain randomization
  // Initialize sum to zero
  dx.zero();
  const atlas::FieldSet zeroFset = util::copyFieldSet(dx[0].fieldSet());

  // Loop over components for the central block
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Randomize covariance
    atlas::FieldSet fset = util::copyFieldSet(zeroFset);
    hybridBlockChain_[jj]->randomize(fset);

    // Apply weight
    // Weight square-root multiplication
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      util::multiplyFieldSet(fset, hybridScalarWeightSqrt_[jj]);
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      util::multiplyFieldSets(fset, hybridFieldWeightSqrt_[jj]);
    }

    // Add component
    util::addFieldSets(dx[0].fieldSet(), fset);
  }

  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(dx[0].fieldSet());

  // ATLAS fieldset to Increment_
  dx[0].synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  // Copy input
  dxo = dxi;

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(dxo[0].fieldSet());

  // Create input FieldSet
  atlas::FieldSet inputFset = util::copyFieldSet(dxo[0].fieldSet());

  // Initialize sum to zero
  util::zeroFieldSet(dxo[0].fieldSet());

  // Loop over B components
  for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
    // Create temporary FieldSet
    atlas::FieldSet fset = util::copyFieldSet(inputFset);

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      util::multiplyFieldSet(fset, hybridScalarWeightSqrt_[jj]);
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      util::multiplyFieldSets(fset, hybridFieldWeightSqrt_[jj]);
    }

    // Apply covariance
    hybridBlockChain_[jj]->multiply(fset);

    // Apply weight
    if (hybridScalarWeightSqrt_[jj] != 1.0) {
      // Scalar weight
      util::multiplyFieldSet(fset, hybridScalarWeightSqrt_[jj]);
    }
    if (!hybridFieldWeightSqrt_[jj].empty()) {
      // File-based weight
      util::multiplyFieldSets(fset, hybridFieldWeightSqrt_[jj]);
    }

    // Add component
    util::addFieldSets(dxo[0].fieldSet(), fset);
  }

  // Apply outer blocks forward
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(dxo[0].fieldSet());

  // ATLAS fieldset to Increment_
  dxo[0].synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
