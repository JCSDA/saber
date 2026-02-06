/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/interface/ModelData.h"
#include "oops/util/ConfigHelpers.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

class SaberParametricBlockChainParameters: public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberParametricBlockChainParameters,
                           ErrorCovarianceParametersBase)

 public:
  // Central and outer blocks
  oops::RequiredParameter<SaberCentralBlockParametersWrapper>
    saberCentralBlockParams{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocksParams{"saber outer blocks", this};

  // Time covariance mode (by default duplicated multivariate)
  // Options: univariate, duplicated multivariate.
  oops::Parameter<std::string> timeCovariance{"time covariance", "multivariate duplicated",
                                              this};

  // Ensemble
  oops::Parameter<bool> iterativeEnsembleLoading{"iterative ensemble loading", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePertOtherGeom{
                        "ensemble pert on other geometry", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleGeom{
                        "ensemble geometry", this};

  // Output ensemble
  oops::OptionalParameter<eckit::LocalConfiguration> outputEnsemble{"output ensemble", this};
};

/// Chain of outer (optional) and not-ensemble central block. Can be used
/// as static error covariance component and as localization for ensemble
/// error covariance.
class SaberParametricBlockChain : public SaberBlockChainBase {
 public:
  /// @brief Standard constructor using MODEL geometry
  template<typename MODEL>
  SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                            const oops::Variables & outerVars,
                            oops::FieldSet4D & fset4dXb,
                            oops::FieldSet4D & fset4dFg,
                            const eckit::Configuration & conf);
  /// @brief Simpler, limited constructor using only generic GeometryData
  SaberParametricBlockChain(const oops::GeometryData & outerGeometryData,
                            const oops::Variables & outerVars,
                            oops::FieldSet4D & fset4dXb,
                            oops::FieldSet4D & fset4dFg,
                            const eckit::Configuration & conf);
  ~SaberParametricBlockChain() = default;

  /// @brief Filter the increment
  void filter(oops::FieldSet4D &) const;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(oops::FieldSet4D &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(oops::FieldSet4D &) const;
  /// @brief Get this B matrix square-root control vector size.
  size_t ctlVecSize() const;
  /// @brief Multiply the control vector by this B matrix square-root.
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const;
  /// @brief Multiply the increment by this B matrix square-root adjoint.
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const;

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const {return outerVariables_;}

 private:
  /// @brief Initialize central block, central function space and central variables.
  ///        Used in constructors.
  std::tuple<oops::Variables, oops::Variables>
      initCentralBlock(const oops::GeometryData & outerGeom,
                       const eckit::Configuration & conf,
                       const SaberBlockParametersBase & saberCentralBlockParams,
                       const oops::FieldSet4D & fset4dXb,
                       const oops::FieldSet4D & fset4dFg);

  /// @brief Run adjoint and square-root tests on central block. Used in constructors.
  void testCentralBlock(const eckit::Configuration & conf,
                        const SaberBlockParametersBase & saberCentralBlockParams,
                        const oops::GeometryData & outerGeom,
                        const oops::Variables & activeVars) const;

  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  bool crossTimeCov_;
  std::unique_ptr<SaberCentralBlockBase> centralBlock_;
  const eckit::mpi::Comm & timeComm_;
  size_t size4D_;
  oops::Variables centralVars_;
  atlas::FunctionSpace centralFunctionSpace_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberParametricBlockChain::SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars),
  timeComm_(fset4dXb.commTime()), size4D_(fset4dXb.size()) {
  oops::Log::trace() << "SaberParametricBlockChain ctor starting" << std::endl;

  // Deserialize parameters and fill configuration with missing values
  SaberParametricBlockChainParameters params;
  params.deserialize(conf);
  eckit::LocalConfiguration fullConf;
  params.serialize(fullConf);

  // Set cross-time covariance flag
  crossTimeCov_ = (params.timeCovariance.value() == "multivariate duplicated");

  // Get central block parameters
  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper
    = params.saberCentralBlockParams;
  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;

  const bool centralDirectCalibration = saberCentralBlockParams.doCalibration();

  // Read ensemble (for non-iterative ensemble loading)
  std::shared_ptr<oops::FieldSets> fsetEns = std::make_shared<oops::FieldSets>(readEnsemble(geom,
                                         outerVars,
                                         fset4dXb.times(), fset4dXb.commTime(), fset4dXb.commEns(),
                                         fullConf));

  // If needed create outer block chain
  if (params.saberOuterBlocksParams.value()) {
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVariables_,
                          fset4dXb, fset4dFg, fullConf,
                          *params.saberOuterBlocksParams.value(),
                          fsetEns, centralDirectCalibration);
  }

  // Set outer geometry data for central block
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geom.generic();


  // Create central block
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  const auto[currentOuterVars, activeVars]
              = initCentralBlock(currentOuterGeom,
                                 fullConf,
                                 saberCentralBlockParams,
                                 fset4dXb,
                                 fset4dFg);

  // Read and add model fields
  centralBlock_->read(geom, currentOuterVars);

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = fullConf.getBool("iterative ensemble loading");

  // Ensemble configuration
  if (saberCentralBlockParams.doCalibration()) {
    // Block calibration
    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

      // Get ensemble size
      size_t nens = getNensFromConfig(fullConf);

      for (size_t ie = 0; ie < nens; ++ie) {
        // Read ensemble member
        oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
        readEnsembleMember(geom, outerVariables_, fullConf, ie, fset);

        // Apply outer blocks inverse (all of them)
        oops::Log::info() << "Info     : Apply outer blocks inverse (all of them)" << std::endl;
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(fset);

        // Use FieldSet in the central block
        oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
        centralBlock_->iterativeCalibrationUpdate(fset);
      }

      // Finalization
      oops::Log::info() << "Info     : Finalization" << std::endl;
      centralBlock_->iterativeCalibrationFinal();
    } else {
      // Direct calibration
      oops::Log::info() << "Info     : Direct calibration" << std::endl;
      centralBlock_->directCalibration(*fsetEns);
    }
  } else if (saberCentralBlockParams.doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    centralBlock_->read();
  }

  if (saberCentralBlockParams.forceWrite.value() || saberCentralBlockParams.doCalibration()) {
    // Write data
    oops::Log::info() << "Info     : Write data" << std::endl;
    centralBlock_->write(geom);
    centralBlock_->write();
  }

  // Write final ensemble
  if (fullConf.has("output ensemble")) {
    // Get output parameters configuration
    const eckit::LocalConfiguration outputEnsembleConf(fullConf, "output ensemble");

    // Check whether geometry grid is similar to the last outer block inner geometry
    const bool useModelWriter = (util::getGridUid(geom.functionSpace())
      == util::getGridUid(currentOuterGeom.functionSpace()));

    // Get ensemble size
    size_t ensembleSize = getNensFromConfig(fullConf);

    // Estimate mean
    oops::FieldSet3D fsetMean(fset4dXb[0].validTime(), geom.getComm());
    if (iterativeEnsembleLoading) {
      for (size_t ie = 0; ie < ensembleSize; ++ie) {
        // Read member
        oops::FieldSet3D fsetMem(fset4dXb[0].validTime(), geom.getComm());
        readEnsembleMember(geom, activeVars, fullConf, ie, fsetMem);

        // Update mean
        if (ie == 0) {
          fsetMean.deepCopy(fsetMem);
        } else {
          fsetMean += fsetMem;
        }
      }

      // Normalize mean
      fsetMean *= 1.0/static_cast<double>(ensembleSize);
    }

    // Write first member only
    const bool firstMemberOnly = outputEnsembleConf.getBool("first member only", false);
    if (firstMemberOnly) {
      ensembleSize = 1;
    }

    for (size_t ie = 0; ie < ensembleSize; ++ie) {
      oops::Log::info() << "Info     : Write member " << ie << std::endl;

      // Increment pointer
      oops::Increment<MODEL> dx(geom, activeVars, fset4dXb[0].validTime());

      // Get ensemble member
      if (iterativeEnsembleLoading) {
        // Read ensemble member
        oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
        readEnsembleMember(geom, activeVars, fullConf, ie, fset);

        // Remove mean
        fset -= fsetMean;

        // Apply outer blocks inverse
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(fset);

        // ATLAS fieldset to Increment_
        dx.fromFieldSet(fset.fieldSet());
      } else {
        // ATLAS fieldset to Increment_
        dx.fromFieldSet((*fsetEns)[ie].fieldSet());
      }

      if (useModelWriter) {
        // Use model writer

        // Set member index
        eckit::LocalConfiguration outputMemberConf(outputEnsembleConf);
        util::setMember(outputMemberConf, ie+1);

        // Write Increment
        dx.write(outputMemberConf);
        oops::Log::test() << "Norm of ensemble member " << ie << ": " << dx.norm() << std::endl;
      } else {
        // Use generic ATLAS writer
        throw eckit::NotImplemented("generic output ensemble write not implemented yet", Here());
      }
    }
  }

  testCentralBlock(fullConf, saberCentralBlockParams, currentOuterGeom, activeVars);

  oops::Log::trace() << "SaberParametricBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
