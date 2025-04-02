/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
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

/// Chain of outer (optional) and not-ensemble central block. Can be used
/// as static error covariance component and as localization for ensemble
/// error covariance.
class SaberParametricBlockChain : public SaberBlockChainBase {
 public:
  /// @brief Standard constructor using MODEL geometry
  template<typename MODEL>
  SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                            const oops::Geometry<MODEL> & dualResGeom,
                            const oops::Variables & outerVars,
                            oops::FieldSet4D & fset4dXb,
                            oops::FieldSet4D & fset4dFg,
                            oops::FieldSets & fsetEns,
                            oops::FieldSets & fsetDualResEns,
                            const eckit::LocalConfiguration & covarConf,
                            const eckit::Configuration & conf);
  /// @brief Simpler, limited constructor using only generic GeometryData
  SaberParametricBlockChain(const oops::GeometryData & outerGeometryData,
                            const oops::Variables & outerVars,
                            oops::FieldSet4D & fset4dXb,
                            oops::FieldSet4D & fset4dFg,
                            const eckit::LocalConfiguration & covarConf,
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
                       const eckit::LocalConfiguration & covarConf,
                       const SaberBlockParametersBase & saberCentralBlockParams,
                       const oops::FieldSet4D & fset4dXb,
                       const oops::FieldSet4D & fset4dFg);

  /// @brief Run adjoint and square-root tests on central block. Used in constructors.
  void testCentralBlock(const eckit::LocalConfiguration & covarConf,
                        const SaberBlockParametersBase & saberCentralBlockParams,
                        const oops::GeometryData & outerGeom,
                        const oops::Variables & activeVars) const;

  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  const bool crossTimeCov_;
  std::unique_ptr<SaberCentralBlockBase> centralBlock_;
  const eckit::mpi::Comm & timeComm_;
  size_t size4D_;
  oops::Variables centralVars_;
  atlas::FunctionSpace centralFunctionSpace_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberParametricBlockChain::SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResGeom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       // TODO(AS): read inside the block so there is no need to pass
                       // as non-const
                       oops::FieldSets & fsetEns,
                       oops::FieldSets & fsetDualResEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf)
  : outerFunctionSpace_(geom.functionSpace()), outerVariables_(outerVars),
  crossTimeCov_(covarConf.getString("time covariance") == "multivariate duplicated"),
  timeComm_(fset4dXb.commTime()), size4D_(fset4dXb.size()) {
  oops::Log::trace() << "SaberParametricBlockChain ctor starting" << std::endl;

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

  // Set outer geometry data for central block
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geom.generic();

  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  const auto[currentOuterVars, activeVars]
              = initCentralBlock(currentOuterGeom,
                                 conf,
                                 covarConf,
                                 saberCentralBlockParams,
                                 fset4dXb,
                                 fset4dFg);

  // Read and add model fields
  centralBlock_->read(geom, currentOuterVars);

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Ensemble configuration
  eckit::LocalConfiguration ensembleConf
         = covarConf.getSubConfiguration("ensemble configuration");
  if (saberCentralBlockParams.doCalibration()) {
    // Block calibration
    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

      // Get ensemble size
      size_t nens = ensembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < nens; ++ie) {
        // Read ensemble member
        oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
        readEnsembleMember(geom, outerVariables_, ensembleConf, ie, fset);

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
      centralBlock_->directCalibration(fsetEns);
    }
  } else if (saberCentralBlockParams.doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    centralBlock_->read();

    if (saberCentralBlockParams.forceWrite.value()) {
      // Write data
      oops::Log::info() << "Info     : Write data" << std::endl;
      centralBlock_->write(geom);
      centralBlock_->write();
    }
  }

  // Dual resolution ensemble
  if (covarConf.has("dual resolution ensemble configuration")) {
    oops::Log::info() << "Info     : Dual resolution setup" << std::endl;

    // Dual resolution setup
    centralBlock_->dualResolutionSetup(dualResGeom.generic());

    // Ensemble configuration
    eckit::LocalConfiguration dualResEnsembleConf
      = covarConf.getSubConfiguration("dual resolution ensemble configuration");

    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

      // Get dual resolution ensemble size
      const size_t dualResNens = dualResEnsembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < dualResNens; ++ie) {
        // Read ensemble member
        oops::FieldSet3D fset(fset4dXb[0].validTime(), dualResGeom.getComm());
        readEnsembleMember(dualResGeom, outerVariables_, dualResEnsembleConf, ie, fset);

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
      centralBlock_->directCalibration(fsetDualResEns);
    }
  }

  // Write calibration data
  if (saberCentralBlockParams.doCalibration()) {
    oops::Log::info() << "Info     : Write calibration data" << std::endl;
    centralBlock_->write(geom);
    centralBlock_->write();
  }

  // Write final ensemble
  if (covarConf.has("output ensemble")) {
    // Get output parameters configuration
    const eckit::LocalConfiguration outputEnsembleConf(covarConf, "output ensemble");

    // Check whether geometry grid is similar to the last outer block inner geometry
    const bool useModelWriter = (util::getGridUid(geom.functionSpace())
      == util::getGridUid(currentOuterGeom.functionSpace()));

    // Get ensemble size
    size_t ensembleSize = ensembleConf.getInt("ensemble size");

    // Estimate mean
    oops::FieldSet3D fsetMean(fset4dXb[0].validTime(), geom.getComm());
    if (iterativeEnsembleLoading) {
      for (size_t ie = 0; ie < ensembleSize; ++ie) {
        // Read member
        oops::FieldSet3D fsetMem(fset4dXb[0].validTime(), geom.getComm());
        readEnsembleMember(geom, activeVars, ensembleConf, ie, fsetMem);

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
        readEnsembleMember(geom, activeVars, ensembleConf, ie, fset);

        // Remove mean
        fset -= fsetMean;

        // Apply outer blocks inverse
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(fset);

        // ATLAS fieldset to Increment_
        dx.fromFieldSet(fset.fieldSet());
      } else {
        // ATLAS fieldset to Increment_
        dx.fromFieldSet(fsetEns[ie].fieldSet());
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

  testCentralBlock(covarConf, saberCentralBlockParams, currentOuterGeom, activeVars);

  oops::Log::trace() << "SaberParametricBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
