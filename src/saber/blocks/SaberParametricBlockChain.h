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

#include "oops/interface/ModelData.h"

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
  template<typename MODEL>
  SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                        const oops::Geometry<MODEL> & dualResGeom,
                        const oops::Variables & outerVars,
                        const atlas::FieldSet & fsetXb,
                        const atlas::FieldSet & fsetFg,
                        const util::DateTime & validTimeOfXbFg,
                        std::vector<atlas::FieldSet> & fsetEns,
                        std::vector<atlas::FieldSet> & dualResolutionFsetEns,
                        const eckit::LocalConfiguration & covarConf,
                        const eckit::Configuration & conf);
  ~SaberParametricBlockChain() = default;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(atlas::FieldSet &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(atlas::FieldSet &) const;

 private:
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  std::unique_ptr<SaberCentralBlockBase> centralBlock_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberParametricBlockChain::SaberParametricBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Geometry<MODEL> & dualResolutionGeom,
                       const oops::Variables & outerVars,
                       const atlas::FieldSet & fsetXb,
                       const atlas::FieldSet & fsetFg,
                       const util::DateTime & validTimeOfXbFg,
                       // TODO(AS): read inside the block so there is no need to pass
                       // as non-const
                       std::vector<atlas::FieldSet> & fsetEns,
                       std::vector<atlas::FieldSet> & dualResolutionFsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const eckit::Configuration & conf) {
  oops::Log::trace() << "SaberParametricBlockChain ctor starting" << std::endl;

  // If needed create outer block chain
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

  // Set outer variables and geometry data for central block
  oops::Variables currentOuterVars = outerBlockChain_ ?
                             outerBlockChain_->innerVars() : outerVars;
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geom.generic();

  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Get active variables
  oops::Variables activeVars = getActiveVars(saberCentralBlockParams, currentOuterVars);
  // Check that active variables are present in variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(currentOuterVars.has(var));
  }

  // Create central block
  centralBlock_ = SaberCentralBlockFactory::create(currentOuterGeom,
                                                   activeVars,
                                                   covarConf,
                                                   saberCentralBlockParams,
                                                   fsetXb,
                                                   fsetFg,
                                                   validTimeOfXbFg,
                                                   geom.timeComm().rank());

  // Read and add model fields
  centralBlock_->read(geom, currentOuterVars, validTimeOfXbFg);

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
        atlas::FieldSet fset;
        readEnsembleMember(geom,
                           outerVars,
                           validTimeOfXbFg,
                           ensembleConf,
                           ie,
                           fset);

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
  }

  // Dual resolution ensemble
  if (covarConf.has("dual resolution ensemble configuration")) {
    oops::Log::info() << "Info     : Dual resolution setup" << std::endl;

    // Dual resolution setup
    centralBlock_->dualResolutionSetup(dualResolutionGeom.generic());

    // Ensemble configuration
    eckit::LocalConfiguration dualResolutionEnsembleConf
      = covarConf.getSubConfiguration("dual resolution ensemble configuration");

    if (iterativeEnsembleLoading) {
      // Iterative calibration
      oops::Log::info() << "Info     : Iterative calibration" << std::endl;

      // Initialization
      centralBlock_->iterativeCalibrationInit();

     // Get dual resolution ensemble size
      size_t dualResolutionNens = dualResolutionEnsembleConf.getInt("ensemble size");

      for (size_t ie = 0; ie < dualResolutionNens; ++ie) {
        // Read ensemble member
        atlas::FieldSet fset;
        readEnsembleMember(dualResolutionGeom,
                           outerVars,
                           validTimeOfXbFg,
                           dualResolutionEnsembleConf,
                           ie,
                           fset);

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
      centralBlock_->directCalibration(dualResolutionFsetEns);
    }
  }

  // Write calibration data
  if (saberCentralBlockParams.doCalibration()) {
    oops::Log::info() << "Info     : Write calibration data" << std::endl;
    centralBlock_->write(geom, currentOuterVars, validTimeOfXbFg);
    centralBlock_->write();
  }

  // Write final ensemble
  if (covarConf.has("output ensemble")) {
    // Write ensemble
    OutputEnsembleParameters<MODEL> outputEnsemble;
    outputEnsemble.deserialize(covarConf.getSubConfiguration("output ensemble"));

    // Check whether geometry grid is similar to the last outer block inner geometry
    const bool useModelWriter = (util::getGridUid(geom.generic().functionSpace())
      == util::getGridUid(currentOuterGeom.functionSpace()));

    // Get ensemble size
    size_t ensembleSize = ensembleConf.getInt("ensemble size");

    // Estimate mean
    atlas::FieldSet mean;
    if (iterativeEnsembleLoading) {
      for (size_t ie = 0; ie < ensembleSize; ++ie) {
        // Read member
        atlas::FieldSet fset;
        readEnsembleMember(geom, activeVars, validTimeOfXbFg, ensembleConf, ie, fset);

        // Update mean
        if (ie == 0) {
          mean = util::copyFieldSet(fset);
        } else {
          util::addFieldSets(mean, fset);
        }
      }

      // Normalize mean
      util::multiplyFieldSet(mean, 1.0/static_cast<double>(ensembleSize));
    }

    // Output parameters
    typename oops::Increment<MODEL>::WriteParameters_ writeParams
      = outputEnsemble.writeParams;
    const bool firstMemberOnly = outputEnsemble.firstMemberOnly.value();

    // Write first member only
    if (firstMemberOnly) {
      ensembleSize = 1;
    }

    for (size_t ie = 0; ie < ensembleSize; ++ie) {
      oops::Log::info() << "Info     : Write member " << ie << std::endl;

      // Increment pointer
      std::unique_ptr<oops::Increment<MODEL>> dx;

      // Get ensemble member
      if (iterativeEnsembleLoading) {
        // Read ensemble member
        dx.reset(new oops::Increment<MODEL>(geom, activeVars, validTimeOfXbFg));
        readEnsembleMember(geom, activeVars, validTimeOfXbFg, ensembleConf, ie, dx->fieldSet());

        // Remove mean
        util::subtractFieldSets(dx->fieldSet(), mean);

        // Apply outer blocks inverse
        if (outerBlockChain_) outerBlockChain_->leftInverseMultiply(dx->fieldSet());
      } else {
        // Copy member
        dx.reset(new oops::Increment<MODEL>(geom, activeVars, validTimeOfXbFg));
        dx->fieldSet() = fsetEns[ie];
      }

      // ATLAS fieldset to Increment_
      dx->synchronizeFields();

      if (useModelWriter) {
        // Use model writer

        // Set member index
        writeParams.setMember(ie+1);

        // Write Increment
        dx->write(writeParams);
        oops::Log::test() << "Norm of ensemble member " << ie << ": " << dx->norm() << std::endl;
      } else {
        // Use generic ATLAS writer
        ABORT("generic output ensemble write not implemented yet");
      }
    }
  }

  // Adjoint test
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    centralBlock_->adjointTest(currentOuterGeom,
                               activeVars,
                               localAdjointTolerance,
                               geom.timeComm().rank());
  }

  oops::Log::trace() << "SaberParametricBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
