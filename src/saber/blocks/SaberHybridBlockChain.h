/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

class CovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters, oops::Parameters)
 public:
  oops::ConfigurationParameter saberBlockChainParams{this};
};

// -----------------------------------------------------------------------------

class WeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WeightParameters, oops::Parameters)
 public:
  // Scalar weight
  oops::Parameter<double> value{"value", 1.0, this};

  // File-base weight
  oops::OptionalParameter<eckit::LocalConfiguration> file{"file", this};
};

// -----------------------------------------------------------------------------

class ComponentParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ComponentParameters, oops::Parameters)
 public:
  // Covariance
  oops::RequiredParameter<CovarianceParameters> covariance{"covariance", this};
  // Weight
  oops::RequiredParameter<WeightParameters> weight{"weight", this};
};

// -----------------------------------------------------------------------------

class SaberHybridBlockChainParameters: public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberHybridBlockChainParameters,
                           ErrorCovarianceParametersBase)
 public:
  // Optional outer blocks
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocksParams{"saber outer blocks", this};
  // Vector of components
  oops::RequiredParameter<std::vector<ComponentParameters>> components{"components", this};
  // Geometry [optional]
  oops::OptionalParameter<eckit::LocalConfiguration> hybridGeometry{"geometry", this};
  // Switch to run components in parallel
  oops::Parameter<bool> runInParallel{"run in parallel", false, this};
  // Switch to run components recursively (for diagnostics)
  oops::Parameter<bool> runComponentsRecursively{"run components recursively", false, this};
  // Resource weighting for each hybrid component.
  oops::OptionalParameter<std::vector<double>> parallelCovarRelativeCPUWeight{
      "parallel covariance relative cpu weight", this};
};

/// Hybrid covariance block chain implementation
template<typename MODEL>
class SaberHybridBlockChain : public SaberBlockChainBase {
 public:
  SaberHybridBlockChain(const oops::Geometry<MODEL> & geom,
                        const oops::Variables & outerVars,
                        oops::FieldSet4D & fset4dXb,
                        oops::FieldSet4D & fset4dFg,
                        const eckit::Configuration & conf);
  ~SaberHybridBlockChain() = default;

  /// @brief Randomize the increment according to this hybrid B matrix.
  void randomize(oops::FieldSet4D &) const override;
  /// @brief Multiply the increment by this hybrid B matrix.
  void multiply(oops::FieldSet4D &) const override;

  /// @brief Control vector size
  size_t ctlVecSize() const override {return 0;}
  /// @brief Square-root multiplication
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const override {}
  /// @brief Adjoint of square-root multiplication
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const override {}

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const override {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const override {return outerVariables_;}

 private:
  /// Function space
  const atlas::FunctionSpace & outerFunctionSpace_;
  /// Variables
  const oops::Variables outerVariables_;

  /// Chain of outer blocks applied to all components of hybrid covariances.
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// Vector of hybrid B components.
  std::vector<std::unique_ptr<SaberBlockChainBase>> hybridBlockChain_;
  /// Vector of scalar weights for hybrid B components.
  std::vector<double> hybridScalarWeightSqrt_;
  /// Vector of field weights for hybrid B components.
  std::vector<oops::FieldSet3D> hybridFieldWeightSqrt_;

  /// Whether to run Hybrid in parallel
  bool parallelHybrid_;
  /// Index of component if running Hybrid in parallel
  size_t myComponent_;
  /// local geometry for parallel Hybrid block (type-erased)
  std::shared_ptr<atlas::FunctionSpace> localHybridFs_, globalHybridFs_;
  std::unique_ptr<oops::Geometry<MODEL>> localHybridGeom_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberHybridBlockChain<MODEL>::SaberHybridBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars),
    parallelHybrid_(false), myComponent_(0) {
  oops::Log::trace() << "SaberHybridBlockChain ctor starting" << std::endl;

  // Deserialize parameters and fill configuration with missing values
  SaberHybridBlockChainParameters params;
  params.deserialize(conf);
  eckit::LocalConfiguration fullConf;
  params.serialize(fullConf);

  // Extract ErrorCovarianceParametersBase from fullConf
  ErrorCovarianceParametersBase paramsBase;
  paramsBase.deserialize(fullConf);

  // Initialize current outer variables
  oops::Variables currentOuterVars(outerVars);

  // Build common (for all hybrid components) outer blocks if they exist
  if (params.saberOuterBlocksParams.value()) {
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVariables_,
                          fset4dXb, fset4dFg, fullConf,
                          *params.saberOuterBlocksParams.value());
    currentOuterVars = outerBlockChain_->innerVars();
  }

  // Hybrid central block
  parallelHybrid_ = params.runInParallel;
  const eckit::mpi::Comm & defaultSpaceComm = geom.getComm();
  const size_t ntasks = defaultSpaceComm.size();
  const size_t nComponents = params.components.value().size();

  std::vector<double> parallelCovRelativeCpuWeight;
  if (params.parallelCovarRelativeCPUWeight.value()) {
    parallelCovRelativeCpuWeight = *params.parallelCovarRelativeCPUWeight.value();
  } else {
    parallelCovRelativeCpuWeight = std::vector<double>(nComponents,
                      1.0 / static_cast<double>(nComponents));
  }

  std::vector<size_t> ntasksPerComponent(nComponents, 0);
  std::vector<size_t> globalTaskOffsetPerComponent(nComponents+1, 0);
  if (parallelHybrid_) {
    oops::Log::info() << "Info     : Creating Hybrid block in parallel" << std::endl;
    // checks
    ASSERT(nComponents == parallelCovRelativeCpuWeight.size());

    // need to check to ensure that the total sum of PEs over components is consistent
    // with the MPI size on the default communicator and that each component
    // has a minimum MPI size of 1.
    for (size_t component = 0; component < nComponents; ++component) {
      ntasksPerComponent[component] =
        std::round(parallelCovRelativeCpuWeight[component] * ntasks);
      ASSERT(ntasksPerComponent[component] > 0);
    }
    int discrepencyPE =
      std::accumulate(ntasksPerComponent.begin(), ntasksPerComponent.end(), 0) - ntasks;

    for (size_t component = 0; component < nComponents && discrepencyPE != 0; ++component) {
      if (discrepencyPE > 0 && ntasksPerComponent[component] >= 2) {
        ntasksPerComponent[component] -= 1;
        discrepencyPE -= 1;
      } else if (discrepencyPE < 0) {
        ntasksPerComponent[component] += 1;
        discrepencyPE += 1;
      }
    }

    ASSERT(std::accumulate(ntasksPerComponent.begin(),
                           ntasksPerComponent.end(), 0) - ntasks == 0);

    for (size_t component = 1; component < nComponents; ++component) {
      globalTaskOffsetPerComponent[component] =
        globalTaskOffsetPerComponent[component-1] + ntasksPerComponent[component-1];
    }
    globalTaskOffsetPerComponent[nComponents] = ntasks;

    const eckit::mpi::Comm & initialDefaultComm = eckit::mpi::comm();
    ASSERT(initialDefaultComm.name() == defaultSpaceComm.name());

    // We split the space communicators only, the time parallelization is untouched
    const size_t myTask = defaultSpaceComm.rank();

    // Set myComponent_  tasksPerComponent
    size_t tasksPerComponent;
    for (size_t component = 0; component < nComponents; ++component) {
      if ((myTask >= globalTaskOffsetPerComponent[component]) &&
          (myTask < globalTaskOffsetPerComponent[component+1])) {
        myComponent_ = component;
        tasksPerComponent = ntasksPerComponent[component];
      }
    }

    oops::Log::info() << "Info     : Creating component " << myComponent_ + 1
                      << "/" << nComponents
                      << " of Hybrid block using " << tasksPerComponent
                      << " MPI tasks." << std::endl;

    // Create communicators for same component, for communications in space
    const auto spaceCommName = ("comm_space_" + std::to_string(myComponent_));
    if (eckit::mpi::hasComm(spaceCommName.c_str())) {
      eckit::mpi::deleteComm(spaceCommName.c_str());
    }
    const auto & localSpaceComm = defaultSpaceComm.split(myComponent_, spaceCommName.c_str());

    // Set up default MPI communicator for atlas
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Create block geometry
    if (params.hybridGeometry.value() == boost::none) {
      throw eckit::UserError("Parallel hybrid block requires geometry key", Here());
    }
    const auto geomConf = *params.hybridGeometry.value();
    // The hybrid Geometry is stored as a class member to ensure it doesn't go
    // out of scope after construction, as it is directly used (not copied) by
    // the hybrid Block Chains.
    localHybridGeom_ =
        std::make_unique<oops::Geometry<MODEL>>(geomConf, localSpaceComm, geom.timeComm());
    localHybridFs_.reset(new atlas::FunctionSpace(localHybridGeom_->functionSpace()));
    globalHybridFs_.reset(new atlas::FunctionSpace(geom.functionSpace()));
    // Copy and redistribute the background and first guess
    oops::FieldSet4D localFset4dXb(fset4dXb.times(), fset4dXb.commTime(),
                                   localSpaceComm, *localHybridFs_, fset4dXb.variables());
    oops::FieldSet4D localFset4dFg(fset4dFg.times(), fset4dFg.commTime(),
                                   localSpaceComm, *localHybridFs_, fset4dFg.variables());

    for (size_t jtime = 0; jtime < fset4dXb.size(); jtime++) {
      util::redistributeToSubcommunicator(fset4dXb[jtime].fieldSet(),
                                          localFset4dXb[jtime].fieldSet(),
                                          defaultSpaceComm,
                                          localSpaceComm,
                                          geom.functionSpace(),
                                          *localHybridFs_);
      util::redistributeToSubcommunicator(fset4dFg[jtime].fieldSet(),
                                          localFset4dFg[jtime].fieldSet(),
                                          defaultSpaceComm,
                                          localSpaceComm,
                                          geom.functionSpace(),
                                          *localHybridFs_);
    }
    defaultSpaceComm.barrier();

    const auto & cmpParams = params.components.value()[myComponent_];

    // Initialize component outer variables
    const oops::Variables cmpOuterVars(currentOuterVars);

    // Set weight
    const auto & weightParams = cmpParams.weight.value();
    // Scalar weight
    hybridScalarWeightSqrt_.push_back(std::sqrt(weightParams.value.value()));
    // File-base weight
    oops::FieldSet3D fsetWeight(localFset4dXb[0].validTime(), localSpaceComm);
    if (weightParams.file.value()) {
      // File-base weight
      readHybridWeight(*localHybridGeom_,
                       cmpOuterVars,
                       localFset4dXb[0].validTime(),
                       *weightParams.file.value(),
                       fsetWeight);
      fsetWeight.sqrt();
    }
    hybridFieldWeightSqrt_.push_back(fsetWeight);

    // Set covariance parameters
    const auto & cmpCovParams = cmpParams.covariance.value();

    // Merge component configuration with full configuration base (order of arguments matters!)
    const eckit::LocalConfiguration cmpMergedConf =
      util::mergeConfigs(cmpCovParams.toConfiguration(), paramsBase.toConfiguration());

    // Add block chain
    hybridBlockChain_.push_back(
        SaberBlockChainFactory<MODEL>::create
         (*localHybridGeom_,
          cmpOuterVars,
          localFset4dXb,
          localFset4dFg,
          cmpMergedConf));

    ASSERT(hybridBlockChain_.size() > 0);

    // Restore previous default MPI communicator for atlas
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());
  } else {
    oops::Log::info() << "Info     : Creating Hybrid block serially" << std::endl;
    // Create block geometry
    const oops::Geometry<MODEL> * hybridGeom = &geom;
    if (params.hybridGeometry.value()) {
      hybridGeom = new oops::Geometry<MODEL>(
        *params.hybridGeometry.value(),
        geom.getComm());
    }
    for (const auto & cmpParams : params.components.value()) {
      // Initialize component outer variables
      const oops::Variables cmpOuterVars(currentOuterVars);

      // Set weight
      const auto & weightParams = cmpParams.weight.value();
      // Scalar weight
      hybridScalarWeightSqrt_.push_back(std::sqrt(weightParams.value));
      // File-base weight
      oops::FieldSet3D fsetWeight(fset4dXb[0].validTime(), geom.getComm());
      if (weightParams.file.value()) {
        // File-base weight
        readHybridWeight(*hybridGeom,
                         cmpOuterVars,
                         fset4dXb[0].validTime(),
                         *weightParams.file.value(),
                         fsetWeight);
        fsetWeight.sqrt();
      }
      hybridFieldWeightSqrt_.push_back(fsetWeight);

      // Set covariance parameters
      const auto & cmpCovParams = cmpParams.covariance.value();

      // Merge component configuration with full configuration base (order of arguments matters!)
      const eckit::LocalConfiguration cmpMergedConf =
        util::mergeConfigs(cmpCovParams.toConfiguration(), paramsBase.toConfiguration());

      // Add block chain
      hybridBlockChain_.push_back
          (SaberBlockChainFactory<MODEL>::create
           (*hybridGeom,
            cmpOuterVars,
            fset4dXb,
            fset4dFg,
            cmpMergedConf));
    }
    ASSERT(hybridBlockChain_.size() > 0);
  }

  oops::Log::trace() << "SaberHybridBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberHybridBlockChain<MODEL>::randomize(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberHybridBlockChain::randomize starting" << std::endl;
  util::Timer timer("SaberHybridBlockChain", "randomize");

  // Initialize FieldSet4D
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    if (outerBlockChain_) {
      fset4d[jtime].init(outerBlockChain_->innerGeometryData().functionSpace(),
        outerBlockChain_->innerVars());
    } else {
      fset4d[jtime].init(outerFunctionSpace_, outerVariables_);
    }
  }
  fset4d.zero();

  if (parallelHybrid_) {
    // Run components of the central block in parallel
    oops::Log::debug() << "Parallel execution of doRandomize in Hybrid" << std::endl;
    oops::Log::debug() << "Running Hybrid component " << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);

    // global communicator and functionSpace
    const auto & defaultSpaceComm = fset4d[0].commGeom();

    // check global communicator is the default one for atlas MPI
    ASSERT(eckit::mpi::comm().name() == defaultSpaceComm.name());

    // subcommunicator within this component
    const auto spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Set up atlas MPI
    eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

    // Create temporary FieldSet on subcommunicator
    oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), localSpaceComm);

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
    defaultSpaceComm.barrier();

    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      // Redistribute to global communicator and sum
       util::gatherAndSumFromSubcommunicator(fset4dCmp[jtime].fieldSet(),
                                             fset4d[jtime].fieldSet(),
                                             localSpaceComm,
                                             defaultSpaceComm,
                                             *localHybridFs_,
                                             *globalHybridFs_);
    }

    // Restore atlas MPI to previous
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());

    fset4d += fset4dCmp;
  } else {
    // Loop over components for the central block
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Randomize covariance
      oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), fset4d[0].commGeom());
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
      fset4d += fset4dCmp;
    }
  }

  if (outerBlockChain_) outerBlockChain_->applyOuterBlocks(fset4d);

  oops::Log::trace() << "SaberHybridBlockChain::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberHybridBlockChain<MODEL>::multiply(oops::FieldSet4D & fset4d) const {
  oops::Log::trace() << "SaberHybridBlockChain::multiply starting" << std::endl;
  util::Timer timer("SaberHybridBlockChain", "multiply");

  // Apply outer blocks adjoint
  if (outerBlockChain_) outerBlockChain_->applyOuterBlocksAD(fset4d);

  // Initialize sum to zero
  oops::FieldSet4D fset4dSum = oops::copyFieldSet4D(fset4d);
  fset4dSum.zero();

  // Loop over B components
  if (parallelHybrid_) {
    oops::Log::debug() << "Parallel execution of Hybrid::multiply, component "
                       << myComponent_ + 1 << std::endl;
    ASSERT(hybridBlockChain_.size() == 1);
    ASSERT(hybridScalarWeightSqrt_.size() == 1);
    ASSERT(hybridFieldWeightSqrt_.size() == 1);

    // Global communicator
    const auto & defaultSpaceComm = fset4d[0].commGeom();
    ASSERT(defaultSpaceComm.name() == eckit::mpi::comm().name());

    // Subcommunicator within component
    const std::string spaceCommName = "comm_space_" + std::to_string(myComponent_);
    const auto & localSpaceComm = eckit::mpi::comm(spaceCommName.c_str());

    // Create temporary FieldSet copy on communicator of this component
    oops::FieldSet4D fset4dCmp(fset4d.times(), fset4d.commTime(), localSpaceComm);
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::redistributeToSubcommunicator(fset4d[jtime].fieldSet(),
                                          fset4dCmp[jtime].fieldSet(),
                                          defaultSpaceComm,
                                          localSpaceComm,
                                          *globalHybridFs_,
                                          *localHybridFs_);
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
    defaultSpaceComm.barrier();

    // Gather and sum data across components
    for (size_t jtime = 0; jtime < fset4dCmp.size(); jtime++) {
      util::gatherAndSumFromSubcommunicator(fset4dCmp[jtime].fieldSet(),
                                            fset4dSum[jtime].fieldSet(),
                                            localSpaceComm,
                                            defaultSpaceComm,
                                            *localHybridFs_,
                                            *globalHybridFs_);
    }

    // Set back default MPI communicator
    eckit::mpi::setCommDefault(defaultSpaceComm.name().c_str());

  } else {
    if (hybridBlockChain_.size() > 1) {
        oops::Log::debug() << "Serial execution of Hybrid::multiply" << std::endl;
    }
    for (size_t jj = 0; jj < hybridBlockChain_.size(); ++jj) {
      // Create temporary FieldSet
      oops::FieldSet4D fset4dCmp = oops::copyFieldSet4D(fset4d);

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

  fset4d.deepCopy(fset4dSum);

  oops::Log::trace() << "SaberHybridBlockChain::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
