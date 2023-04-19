/*
 * (C) Copyright 2021 UCAR
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
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/SaberBlockChain.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef oops::State<MODEL>                                   State_;
  typedef typename oops::Increment<MODEL>::WriteParameters_    WriteParameters_;
  typedef typename boost::ptr_vector<SaberBlockChain>          SaberBlockChainVec_;
  typedef typename SaberBlockChainVec_::const_iterator         chainIcst_;

 public:
  typedef ErrorCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const Parameters_ &,
                  const State_ &, const State_ &);
  virtual ~ErrorCovariance();

  // Required by iterative inverse
  void multiply(const Increment_ & dxi, Increment_ & dxo) const {this->doMultiply(dxi, dxo);}

 private:
  ErrorCovariance(const ErrorCovariance&);
  ErrorCovariance& operator=(const ErrorCovariance&);

  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<SaberBlockChain> singleBlockChain_;
  SaberBlockChainVec_ hybridBlockChain_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const Parameters_ & params,
                                        const State_ & xb,
                                        const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, params, xb, fg), singleBlockChain_(),
    hybridBlockChain_()
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;

  // Local copy of background and first guess
  State_ xbLocal(xb);
  State_ fgLocal(fg);

  // Extend backgroud and first guess with extra fields
  // TODO(Benjamin, Marek, Mayeul, ?)

  // Initialize geometry data
  std::vector<std::reference_wrapper<const oops::GeometryData>> outerGeometryData;
  outerGeometryData.push_back(geom.generic());

  // Intialize outer variables
  oops::Variables outerVars(incVars);

  // Initialize single blockchain
  Increment_ dx(geom, incVars, xb.validTime());
  singleBlockChain_.reset(new SaberBlockChain(incVars, dx.fieldSet()));

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = params.iterativeEnsembleLoading.value();

  // Initialize ensembles as vector of FieldSets
  std::vector<atlas::FieldSet> fsetEns;

  // Read ensemble
  eckit::LocalConfiguration ensembleConf = readEnsemble(geom,
                                                        incVars,
                                                        xbLocal,
                                                        fgLocal,
                                                        params.toConfiguration(),
                                                        iterativeEnsembleLoading,
                                                        fsetEns);

  // Create covariance configuration
  eckit::LocalConfiguration covarConf;
  covarConf.set("ensemble configuration", ensembleConf);
  covarConf.set("adjoint test", params.adjointTest.value());
  covarConf.set("adjoint tolerance", params.adjointTolerance.value());
  covarConf.set("inverse test", params.inverseTest.value());
  covarConf.set("inverse tolerance", params.inverseTolerance.value());
  covarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());

  // Build outer blocks successively
  const auto & saberOuterBlocksParams = params.saberOuterBlocksParams.value();
  if (saberOuterBlocksParams != boost::none) {
    buildOuterBlocks(geom,
                     outerGeometryData,
                     outerVars,
                     xbLocal,
                     fgLocal,
                     fsetEns,
                     covarConf,
                     *saberOuterBlocksParams,
                     *singleBlockChain_);
  }

  // Dual resolution ensemble
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
    oops::State<MODEL> xbDualResolution(*dualResolutionGeom, xb);
    oops::State<MODEL> fgDualResolution(*dualResolutionGeom, fg);

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

  // Build central block
  const SaberBlockParametersBase & saberCentralBlockParams =
    params.saberCentralBlockParams.value().saberCentralBlockParameters;
  if (saberCentralBlockParams.saberBlockName.value() == "Hybrid") {
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
      // Initialize component geometry data
      std::vector<std::reference_wrapper<const oops::GeometryData>> cmpOuterGeometryData;
      cmpOuterGeometryData.push_back(outerGeometryData.back().get());

      // Intialize component outer variables
      oops::Variables cmpOuterVars(outerVars);

      // Initialize ensembles as vector of FieldSets
      std::vector<atlas::FieldSet> cmpFsetEns;

      // Initialize hybrid blockchain component
      hybridBlockChain_.push_back(new SaberBlockChain(cmpOuterVars,
        singleBlockChain_->centralFieldSet()));

      // Set weight
      eckit::LocalConfiguration weightConf = cmp.getSubConfiguration("weight");
      if (weightConf.has("value")) {
        // Scalar weight
        hybridBlockChain_.back().setWeight(weightConf.getDouble("value"));
      } else if (weightConf.has("file")) {
        // File-base weight
        atlas::FieldSet fset;
        readHybridWeight(*hybridGeom,
                         outerVars,
                         xbLocal.validTime(),
                         weightConf.getSubConfiguration("file"),
                         fset);
        hybridBlockChain_.back().setWeight(fset);
      } else {
        ABORT("missing hybrid weight");
      }

      // Set covariance
      eckit::LocalConfiguration cmpConf = cmp.getSubConfiguration("covariance");

      // Read ensemble
      eckit::LocalConfiguration cmpEnsembleConf
         = readEnsemble(*hybridGeom,
                        cmpOuterVars,
                        xbLocal,
                        fgLocal,
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

      // Build outer blocks successively
      if (cmpConf.has("saber outer blocks")) {
        std::vector<SaberOuterBlockParametersWrapper> cmpOuterBlocksParams;
        for (const auto & cmpOuterBlockConf : cmpConf.getSubConfigurations("saber outer blocks")) {
          SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
          cmpOuterBlockParamsWrapper.deserialize(cmpOuterBlockConf);
          cmpOuterBlocksParams.push_back(cmpOuterBlockParamsWrapper);
        }
        buildOuterBlocks(*hybridGeom,
                         cmpOuterGeometryData,
                         cmpOuterVars,
                         xbLocal,
                         fgLocal,
                         cmpFsetEns,
                         cmpCovarConf,
                         cmpOuterBlocksParams,
                         hybridBlockChain_.back());
      }

      // Central block
      SaberCentralBlockParametersWrapper cmpCentralBlockParamsWrapper;
      cmpCentralBlockParamsWrapper.deserialize(cmpConf.getSubConfiguration("saber central block"));
      buildCentralBlock(*hybridGeom,
                        *dualResolutionGeom,
                        cmpOuterGeometryData.back().get(),
                        cmpOuterVars,
                        xbLocal,
                        fgLocal,
                        cmpFsetEns,
                        dualResolutionFsetEns,
                        cmpCovarConf,
                        cmpCentralBlockParamsWrapper,
                        hybridBlockChain_.back());
    }
  } else {
    // Single central block
    buildCentralBlock(geom,
                      *dualResolutionGeom,
                      outerGeometryData.back().get(),
                      outerVars,
                      xbLocal,
                      fgLocal,
                      fsetEns,
                      dualResolutionFsetEns,
                      covarConf,
                      params.saberCentralBlockParams.value(),
                      *singleBlockChain_);
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
void ErrorCovariance<MODEL>::doRandomize(Increment_ & dx) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  // SABER block chain randomization
  if (hybridBlockChain_.size() > 0) {
    // Hybrid central block

    // Initialize sum to zero
    util::zeroFieldSet(dx.fieldSet());

    // Loop over components for the central block
    for (chainIcst_ it = hybridBlockChain_.begin(); it != hybridBlockChain_.end(); ++it) {
      // Create temporary FieldSet
      atlas::FieldSet fset = util::copyFieldSet(dx.fieldSet());

      // Randomize covariance
      it->randomize(fset);

      // Add component
      util::addFieldSets(dx.fieldSet(), fset);
    }

    // Apply outer blocks forward
    singleBlockChain_->applyOuterBlocks(dx.fieldSet());
  } else {
    // Single central block
    singleBlockChain_->randomize(dx.fieldSet());
  }

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment_ & dxi,
                                        Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  // Copy input
  dxo = dxi;

  // SABER block chain multiplication
  if (hybridBlockChain_.size() > 0) {
    // Hybrid central block

    // Apply outer blocks adjoint
    singleBlockChain_->applyOuterBlocksAD(dxo.fieldSet());

    // Create input FieldSet
    atlas::FieldSet inputFset = util::copyFieldSet(dxo.fieldSet());

    // Initialize sum to zero
    util::zeroFieldSet(dxo.fieldSet());

    // Loop over components for the central block
    for (chainIcst_ it = hybridBlockChain_.begin(); it != hybridBlockChain_.end(); ++it) {
      // Create temporary FieldSet
      atlas::FieldSet fset = util::copyFieldSet(inputFset);

      // Apply covariance
      it->multiply(fset);

      // Add component
      util::addFieldSets(dxo.fieldSet(), fset);
    }

    // Apply outer blocks forward
    singleBlockChain_->applyOuterBlocks(dxo.fieldSet());
  } else {
    // Single central block
    singleBlockChain_->multiply(dxo.fieldSet());
  }

  // ATLAS fieldset to Increment_
  dxo.synchronizeFields();

  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                               Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Iterative inverse
  oops::IdentityMatrix<Increment_> Id;
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
