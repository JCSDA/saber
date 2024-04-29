/*
 * (C) Copyright 2021-2023 UCAR
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
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

inline std::string parametricIfNotEnsemble(const std::string & blockName) {
  if (blockName == "Ensemble")
    return blockName;
  else if (blockName == "gsi hybrid covariance")
    return "GSI";
  else
    return "Parametric";
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class ErrorCovariance : public oops::ModelSpaceCovarianceBase<MODEL>,
                        public util::Printable,
                        private util::ObjectCounter<ErrorCovariance<MODEL>> {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef oops::Increment4D<MODEL>                             Increment4D_;
  typedef oops::State4D<MODEL>                                 State4D_;

 public:
  static const std::string classname() {return "saber::ErrorCovariance";}

  ErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const eckit::Configuration &,
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
  std::vector<oops::FieldSet3D> hybridFieldWeightSqrt_;

  /// Whether to run Hybrid in parallel
  bool parallelHybrid_;
  /// Index of component if running Hybrid in parallel
  size_t myComponent_;  // This is not strictly necessary
  /// local geometry just out of parallel Hybrid block
  std::shared_ptr<Geometry_> localHybridGeom_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance<MODEL>::ErrorCovariance(const Geometry_ & geom,
                                        const oops::Variables & incVars,
                                        const eckit::Configuration & config,
                                        const State4D_ & xb,
                                        const State4D_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(geom, config, xb, fg),
    parallelHybrid_(false),
    myComponent_(-1)
{
  oops::Log::trace() << "ErrorCovariance::ErrorCovariance starting" << std::endl;
  ErrorCovarianceParameters<MODEL> params;
  params.deserialize(config);

  // Local copy of background and first guess that can undergo interpolation
  const oops::FieldSet4D fset4dXbTmp(xb);
  const oops::FieldSet4D fset4dFgTmp(fg);

  oops::FieldSet4D fset4dXb = oops::copyFieldSet4D(fset4dXbTmp);
  oops::FieldSet4D fset4dFg = oops::copyFieldSet4D(fset4dFgTmp);

  // Extend background and first guess with geometry fields
  // TODO(Benjamin, Marek, Mayeul, ?)

  // Initialize outer variables
  const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
  oops::Variables outerVars(incVars);
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    outerVars.addMetaData(outerVars[i], "levels", vlevs[i]);
  }

  // Create covariance configuration
  eckit::LocalConfiguration covarConf;
  covarConf.set("adjoint test", params.adjointTest.value());
  covarConf.set("adjoint tolerance", params.adjointTolerance.value());
  covarConf.set("inverse test", params.inverseTest.value());
  covarConf.set("inverse tolerance", params.inverseTolerance.value());
  covarConf.set("square-root test", params.sqrtTest.value());
  covarConf.set("square-root tolerance", params.sqrtTolerance.value());
  covarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());
  covarConf.set("time covariance", params.timeCovariance.value());

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = params.iterativeEnsembleLoading.value();

  // Read ensemble (for non-iterative ensemble loading)
  eckit::LocalConfiguration ensembleConf;
  oops::FieldSets fsetEns = readEnsemble(geom,
                                         outerVars,
                                         xb,
                                         fg,
                                         params.toConfiguration(),
                                         iterativeEnsembleLoading,
                                         ensembleConf);

  covarConf.set("ensemble configuration", ensembleConf);
  // Read dual resolution ensemble if needed
  const auto & dualResParams = params.dualResParams.value();
  const Geometry_ * dualResGeom = &geom;
  std::unique_ptr<oops::FieldSets> fsetDualResEns;
  if (dualResParams != boost::none) {
    const auto & dualResGeomConf = dualResParams->geometry.value();
    if (dualResGeomConf != boost::none) {
      // Create dualRes geometry
      typename Geometry_::Parameters_ dualResGeomParams;
      dualResGeomParams.deserialize(*dualResGeomConf);
      dualResGeom = new Geometry_(dualResGeomParams, geom.getComm());
    }
    // Background and first guess at dual resolution geometry
    const State4D_ xbDualRes(*dualResGeom, xb);
    const State4D_ fgDualRes(*dualResGeom, fg);
    // Read dual resolution ensemble
    eckit::LocalConfiguration dualResEnsembleConf;
    fsetDualResEns = std::make_unique<oops::FieldSets>(readEnsemble(*dualResGeom,
                     outerVars,
                     xbDualRes,
                     fgDualRes,
                     dualResParams->toConfiguration(),
                     iterativeEnsembleLoading,
                     dualResEnsembleConf));
    // Add dual resolution ensemble configuration
    covarConf.set("dual resolution ensemble configuration", dualResEnsembleConf);
  }
  if (!fsetDualResEns) {
    std::vector<util::DateTime> dates;
    std::vector<int> ensmems;
    fsetDualResEns = std::make_unique<oops::FieldSets>(dates,
                                      xb.commTime(), ensmems, xb.commEns());
  }

  // Add ensemble output
  const auto & outputEnsemble = params.outputEnsemble.value();
  if (outputEnsemble != boost::none) {
    covarConf.set("output ensemble", *outputEnsemble);
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
                       fset4dXb,
                       fset4dFg,
                       fsetEns,
                       covarConf,
                       *saberOuterBlocksParams);
      outerVars = outerBlockChain_->innerVars();
    }

    // Hybrid central block
    eckit::LocalConfiguration hybridConf = saberCentralBlockParams.toConfiguration();

    parallelHybrid_ = hybridConf.getBool("run in parallel");

    const size_t nComponents = hybridConf.getSubConfigurations("components").size();
    const eckit::mpi::Comm & globalSpaceComm = geom.getComm();
    const size_t ntasks = globalSpaceComm.size();

    if (parallelHybrid_ && ntasks % nComponents != 0) {
      oops::Log::warning() << "Warning  : Number of MPI tasks not divisible "
                              "by number of Hybrid block components, running serially."
                           << std::endl;
      parallelHybrid_ = false;
    }

    if (parallelHybrid_) {
      oops::Log::info() << "Info     : Creating Hybrid block in parallel" << std::endl;

      if (dualResParams != boost::none) {
        throw eckit::NotImplemented("Parallel Hybrid not compatible "
                                    "with dual resolution ensemble yet",
                                    Here());
      }

      const eckit::mpi::Comm & initialDefaultComm = eckit::mpi::comm();
      ASSERT(initialDefaultComm.name() == globalSpaceComm.name());

      // We split the space communicators only, the time parallelization is untouched
      const size_t myTask = globalSpaceComm.rank();
      const size_t tasksPerComponent = ntasks / nComponents;
      myComponent_ = myTask / tasksPerComponent;

      oops::Log::info() << "Info     : Creating component " << myComponent_ + 1
                        << "/" << nComponents
                        << " of Hybrid block using " << tasksPerComponent
                        << " MPI tasks." << std::endl;

      // Create communicators for same component, for communications in space
      const auto spaceCommName = ("comm_space_" + std::to_string(myComponent_)).c_str();
      if (eckit::mpi::hasComm(spaceCommName)) {
        eckit::mpi::deleteComm(spaceCommName);
      }
      const auto & localSpaceComm = globalSpaceComm.split(myComponent_, spaceCommName);

      // Create block geometry (needed for ensemble reading and local geometries)
      if (!hybridConf.has("geometry")) {
        throw eckit::UserError("Parallel hybrid block requires geometry key", Here());
      }
      const auto geomConf = hybridConf.getSubConfiguration("geometry");
      // The hybrid Geometry is stored as a class member to ensure it doesn't go
      // out of scope after construction, as it is directly used (not copied) by
      // the hybrid Block Chains.
      localHybridGeom_.reset(new Geometry_(geomConf, localSpaceComm, geom.timeComm()));

      // Copy and redistribute the background and first guess
      State4D_ localXb(*localHybridGeom_, xb.variables(), xb.times(), xb.commTime());
      State4D_ localFg(*localHybridGeom_, fg.variables(), fg.times(), fg.commTime());

      for (size_t jtime = 0; jtime < xb.size(); jtime++) {
        util::redistributeToSubcommunicator(xb[jtime].fieldSet().fieldSet(),
                                            localXb[jtime].fieldSet().fieldSet(),
                                            globalSpaceComm,
                                            localSpaceComm,
                                            geom.functionSpace(),
                                            localHybridGeom_->functionSpace());
        util::redistributeToSubcommunicator(fg[jtime].fieldSet().fieldSet(),
                                            localFg[jtime].fieldSet().fieldSet(),
                                            globalSpaceComm,
                                            localSpaceComm,
                                            geom.functionSpace(),
                                            localHybridGeom_->functionSpace());
      }
      globalSpaceComm.barrier();

      // Set up default MPI communicator for atlas
      eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

      const oops::FieldSet4D localFset4dXbTmp(localXb);
      const oops::FieldSet4D localFset4dFgTmp(localFg);

      oops::FieldSet4D localFset4dXb = oops::copyFieldSet4D(localFset4dXbTmp);
      oops::FieldSet4D localFset4dFg = oops::copyFieldSet4D(localFset4dFgTmp);

      const auto cmp = hybridConf.getSubConfigurations("components")[myComponent_];

      // Initialize component outer variables
      const oops::Variables cmpOuterVars(outerVars);

      // Set weight
      eckit::LocalConfiguration weightConf = cmp.getSubConfiguration("weight");
      // Scalar weight
      hybridScalarWeightSqrt_.push_back(std::sqrt(weightConf.getDouble("value", 1.0)));

      // File-base weight
      oops::FieldSet3D fsetWeight(localXb[0].validTime(), localSpaceComm);
      if (weightConf.has("file")) {
        // File-base weight
        readHybridWeight(*localHybridGeom_,
                         outerVars,
                         localXb[0].validTime(),
                         weightConf.getSubConfiguration("file"),
                         fsetWeight);
        fsetWeight.sqrt();
      }
      hybridFieldWeightSqrt_.push_back(fsetWeight);

      // Set covariance
      eckit::LocalConfiguration cmpConf = cmp.getSubConfiguration("covariance");

      // Read ensemble
      eckit::LocalConfiguration cmpEnsembleConf;
      oops::FieldSets localFset4dCmpEns
         = readEnsemble(*localHybridGeom_,
                        cmpOuterVars,
                        localXb,
                        localFg,
                        cmpConf,
                        params.iterativeEnsembleLoading.value(),
                        cmpEnsembleConf);

      // Create internal configuration
      eckit::LocalConfiguration cmpCovarConf;
      cmpCovarConf.set("ensemble configuration", cmpEnsembleConf);
      cmpCovarConf.set("adjoint test", params.adjointTest.value());
      cmpCovarConf.set("adjoint tolerance", params.adjointTolerance.value());
      cmpCovarConf.set("inverse test", params.inverseTest.value());
      cmpCovarConf.set("inverse tolerance", params.inverseTolerance.value());
      cmpCovarConf.set("square-root test", params.sqrtTest.value());
      cmpCovarConf.set("square-root tolerance", params.sqrtTolerance.value());
      cmpCovarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());
      cmpCovarConf.set("time covariance", params.timeCovariance.value());

      SaberCentralBlockParametersWrapper cmpCentralBlockParamsWrapper;
      cmpCentralBlockParamsWrapper.deserialize(cmpConf.getSubConfiguration("saber central block"));
      const auto & centralBlockParams =
                   cmpCentralBlockParamsWrapper.saberCentralBlockParameters.value();

      hybridBlockChain_.push_back(
        SaberBlockChainFactory<MODEL>::create
         (parametricIfNotEnsemble(centralBlockParams.saberBlockName.value()),
          *localHybridGeom_,
          *dualResGeom,
          cmpOuterVars,
          localFset4dXb,
          localFset4dFg,
          localFset4dCmpEns,
          *fsetDualResEns,
          cmpCovarConf,
          cmpConf));

      ASSERT(hybridBlockChain_.size() > 0);

      // Restore previous default MPI communicator for atlas
      eckit::mpi::setCommDefault(globalSpaceComm.name().c_str());
    } else {
      oops::Log::info() << "Info     : Creating Hybrid block serially" << std::endl;
      // Create block geometry (needed for ensemble reading)
      const Geometry_ * hybridGeom = &geom;
      if (hybridConf.has("geometry")) {
        hybridGeom = new Geometry_(hybridConf.getSubConfiguration("geometry"),
          geom.getComm());
      }
      for (const auto & cmp : hybridConf.getSubConfigurations("components")) {
        // Initialize component outer variables
        const oops::Variables cmpOuterVars(outerVars);

        // Set weight
        eckit::LocalConfiguration weightConf = cmp.getSubConfiguration("weight");
        // Scalar weight
        hybridScalarWeightSqrt_.push_back(std::sqrt(weightConf.getDouble("value", 1.0)));
        // File-base weight
        oops::FieldSet3D fsetWeight(xb[0].validTime(), geom.getComm());
        if (weightConf.has("file")) {
          // File-base weight
          readHybridWeight(*hybridGeom,
                           outerVars,
                           xb[0].validTime(),
                           weightConf.getSubConfiguration("file"),
                           fsetWeight);
          fsetWeight.sqrt();
        }
        hybridFieldWeightSqrt_.push_back(fsetWeight);

        // Set covariance
        eckit::LocalConfiguration cmpConf = cmp.getSubConfiguration("covariance");

        // Read ensemble
        eckit::LocalConfiguration cmpEnsembleConf;
        oops::FieldSets fset4dCmpEns
           = readEnsemble(*hybridGeom,
                          cmpOuterVars,
                          xb,
                          fg,
                          cmpConf,
                          params.iterativeEnsembleLoading.value(),
                          cmpEnsembleConf);

        // Create internal configuration
        eckit::LocalConfiguration cmpCovarConf;
        cmpCovarConf.set("ensemble configuration", cmpEnsembleConf);
        cmpCovarConf.set("adjoint test", params.adjointTest.value());
        cmpCovarConf.set("adjoint tolerance", params.adjointTolerance.value());
        cmpCovarConf.set("inverse test", params.inverseTest.value());
        cmpCovarConf.set("inverse tolerance", params.inverseTolerance.value());
        cmpCovarConf.set("square-root test", params.sqrtTest.value());
        cmpCovarConf.set("square-root tolerance", params.sqrtTolerance.value());
        cmpCovarConf.set("iterative ensemble loading", params.iterativeEnsembleLoading.value());
        cmpCovarConf.set("time covariance", params.timeCovariance.value());

        SaberCentralBlockParametersWrapper cmpCentralBlockParamsWrapper;
        cmpCentralBlockParamsWrapper.deserialize(
                    cmpConf.getSubConfiguration("saber central block"));
        const auto & centralBlockParams =
                     cmpCentralBlockParamsWrapper.saberCentralBlockParameters.value();

        hybridBlockChain_.push_back
          (SaberBlockChainFactory<MODEL>::create
           (parametricIfNotEnsemble(centralBlockParams.saberBlockName.value()),
            *hybridGeom,
            *dualResGeom,
            cmpOuterVars,
            fset4dXb,
            fset4dFg,
            fset4dCmpEns,
            *fsetDualResEns,
            cmpCovarConf,
            cmpConf));
      }
      ASSERT(hybridBlockChain_.size() > 0);
    }
  } else {
    // Non-hybrid covariance: single block chain
    hybridBlockChain_.push_back
      (SaberBlockChainFactory<MODEL>::create
       (parametricIfNotEnsemble(saberCentralBlockParams.saberBlockName.value()),
        geom,
        *dualResGeom,
        outerVars,
        fset4dXb,
        fset4dFg,
        fsetEns,
        *fsetDualResEns,
        covarConf,
        params.toConfiguration()));

    // Set weights
    hybridScalarWeightSqrt_.push_back(1.0);
    // File-base weight
    oops::FieldSet3D fsetWeight(xb[0].validTime(), geom.getComm());
    hybridFieldWeightSqrt_.push_back(fsetWeight);
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

  // Create FieldSet4D, set to zero
  oops::FieldSet4D fset4dSum(dx.times(), dx.commTime(), dx.geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum.size(); ++jtime) {
    fset4dSum[jtime].init(hybridBlockChain_[0]->outerFunctionSpace(),
                          hybridBlockChain_[0]->outerVariables(),
                          0.0);
  }

  if (parallelHybrid_) {
    // Run components of the central block in parallel
    oops::Log::debug() << "Parallel execution of doRandomize in Hybrid" << std::endl;
    oops::Log::debug() << "Running Hybrid component " << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);

    // global communicator and functionSpace
    const auto & globalSpaceComm = dx.geometry().getComm();
    const auto & globalFunctionSpace = dx.geometry().functionSpace();

    // check global communicator is the default one for atlas MPI
    ASSERT(eckit::mpi::comm().name() == globalSpaceComm.name());

    // subcommunicator within this component
    const auto spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Set up atlas MPI
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Create temporary FieldSet on subcommunicator
    oops::FieldSet4D fset4dCmp(dx.times(), dx.commTime(), localSpaceComm);

    hybridBlockChain_[0]->randomize(fset4dCmp);

    // Weight square-root multiplication
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Add components
    globalSpaceComm.barrier();

    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      // Redistribute to global communicator and sum
       util::gatherAndSumFromSubcommunicator(fset4dCmp[jtime].fieldSet(),
                                             fset4dSum[jtime].fieldSet(),
                                             localSpaceComm,
                                             globalSpaceComm,
                                             localHybridGeom_->functionSpace(),
                                             globalFunctionSpace);
    }

    // Restore atlas MPI to previous
    eckit::mpi::setCommDefault(globalSpaceComm.name().c_str());

    fset4dSum += fset4dCmp;
  } else {
    // Loop over components for the central block
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Randomize covariance
      oops::FieldSet4D fset4dCmp(dx.times(), dx.commTime(), dx.geometry().getComm());
      hybridBlockChain_[jj]->randomize(fset4dCmp);

      // Weight square-root multiplication
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Add component
      fset4dSum += fset4dCmp;
    }
  }

  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dx.size(); ++jtime) {
    dx[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

  oops::Log::trace() << "ErrorCovariance<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance<MODEL>::doMultiply(const Increment4D_ & dxi,
                                        Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  // Copy input
  dxo = dxi;
  oops::FieldSet4D fset4dInit(dxi);

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(fset4dInit);

  // Initialize sum to zero
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4dInit);
  fset4dSum.zero();

  // Loop over B components
  if (parallelHybrid_) {
    oops::Log::debug() << "Parallel execution of Hybrid::multiply, component "
                       << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);
    ASSERT(hybridScalarWeightSqrt_.size() == 1);
    ASSERT(hybridFieldWeightSqrt_.size() == 1);

    // Global communicator
    const auto & globalSpaceComm = dxi.geometry().getComm();
    const auto & globalFunctionSpace = dxi.geometry().functionSpace();
    ASSERT(globalSpaceComm.name() == eckit::mpi::comm().name());

    // Subcommunicator within component
    const std::string spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Create temporary FieldSet copy on communicator of this component
    oops::FieldSet4D fset4dCmp(fset4dInit.times(), fset4dInit.commTime(), localSpaceComm);
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::redistributeToSubcommunicator(fset4dInit[jtime].fieldSet(),
                                          fset4dCmp[jtime].fieldSet(),
                                          globalSpaceComm,
                                          localSpaceComm,
                                          globalFunctionSpace,
                                          localHybridGeom_->functionSpace());
    }

    // Set up atlas MPI
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Apply weight
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Apply covariance
    hybridBlockChain_[0]->multiply(fset4dCmp);

    // Apply weight
    if (hybridScalarWeightSqrt_[0] != 1.0) {
      // Scalar weight
      fset4dCmp *= hybridScalarWeightSqrt_[0];
    }
    if (!hybridFieldWeightSqrt_[0].empty()) {
      // File-based weight
      fset4dCmp *= hybridFieldWeightSqrt_[0];
    }

    // Wait for all components to have finished multiplying
    globalSpaceComm.barrier();

    // Gather and sum data across components
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::gatherAndSumFromSubcommunicator(fset4dCmp[jtime].fieldSet(),
                                            fset4dSum[jtime].fieldSet(),
                                            localSpaceComm,
                                            globalSpaceComm,
                                            localHybridGeom_->functionSpace(),
                                            globalFunctionSpace);
    }

    // Set back default MPI communicator
    eckit::mpi::setCommDefault(globalSpaceComm.name().c_str());

  } else {
    if (hybridBlockChain_.size() > 1) {
        oops::Log::debug() << "Serial execution of Hybrid::multiply" << std::endl;
    }
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Create temporary FieldSet
      oops::FieldSet4D fset4dCmp = oops::copyFieldSet4D(fset4dInit);

      // Apply weight
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Apply covariance
      hybridBlockChain_[jj]->multiply(fset4dCmp);

      // Apply weight
      if (hybridScalarWeightSqrt_[jj] != 1.0) {
        // Scalar weight
        fset4dCmp *= hybridScalarWeightSqrt_[jj];
      }
      if (!hybridFieldWeightSqrt_[jj].empty()) {
        // File-based weight
        fset4dCmp *= hybridFieldWeightSqrt_[jj];
      }

      // Add component
      fset4dSum += fset4dCmp;
    }
  }

  // Apply outer blocks forward
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4dSum);

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dxo.size(); ++jtime) {
    dxo[jtime].fromFieldSet(fset4dSum[jtime].fieldSet());
  }

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
