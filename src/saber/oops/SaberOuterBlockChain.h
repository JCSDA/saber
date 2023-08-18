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

#include "oops/base/Geometry.h"
#include "oops/interface/ModelData.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {

/// Chain of outer saber blocks (no central block). Can be used as the common
/// outer blocks for the hybrid covariance, outer blocks for static and ensemble
/// covariances, ensemble transform for the ensemble covariance.
class SaberOuterBlockChain {
 public:
  template<typename MODEL>
  SaberOuterBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       const atlas::FieldSet & fsetXb,
                       const atlas::FieldSet & fsetFg,
                       const util::DateTime & validTimeOfXbFg,
                       std::vector<atlas::FieldSet> & fsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const std::vector<saber::SaberOuterBlockParametersWrapper> & params);
  ~SaberOuterBlockChain() = default;

  // Accessors
  const std::vector<std::unique_ptr<SaberOuterBlockBase>> & outerBlocks() const
    {return outerBlocks_;}
  // TODO(AS): remove non-const accessor (currently used to add ens transform blocks
  // to the outer blocks in SaberEnsembleBlockChain)
  std::vector<std::unique_ptr<SaberOuterBlockBase>> & outerBlocks() {return outerBlocks_;}

  /// @brief Returns inner-most variables.
  const oops::Variables & innerVars() const {
    return outerBlocks_.back()->innerVars();
  }
  /// @brief Returns inner-most geometry data.
  const oops::GeometryData & innerGeometryData() const {
    return outerBlocks_.back()->innerGeometryData();
  }

  /// @brief Forward multiplication by all outer blocks.
  void applyOuterBlocks(atlas::FieldSet & fset) const {
    for (auto it = outerBlocks_.rbegin(); it != outerBlocks_.rend(); ++it) {
      it->get()->multiply(fset);
    }
  }

  /// @brief Adjoint multiplication by all outer blocks.
  void applyOuterBlocksAD(atlas::FieldSet & fset) const {
    // Outer blocks adjoint multiplication
    for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
      it->get()->multiplyAD(fset);
    }
  }

  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the ones that haven't implemented inverse yet.
  void leftInverseMultiply(atlas::FieldSet & fset) const {
    for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
      if (it->get()->skipInverse()) {
        oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                        << it->get()->blockName() << std::endl;
      } else {
        it->get()->leftInverseMultiply(fset);
      }
    }
  }

 private:
  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the last one and the ones that haven't implemented inverse yet.
  void leftInverseMultiplyExceptLast(atlas::FieldSet & fset) const {
    // Outer blocks left inverse multiplication
    for (auto it = outerBlocks_.begin(); it != std::prev(outerBlocks_.end()); ++it) {
      if (it->get()->skipInverse()) {
        oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                        << it->get()->blockName() << std::endl;
      } else {
        it->get()->leftInverseMultiply(fset);
      }
    }
  }

  /// @brief Vector of all outer blocks.
  std::vector<std::unique_ptr<SaberOuterBlockBase>> outerBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberOuterBlockChain::SaberOuterBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       const atlas::FieldSet & fsetXb,
                       const atlas::FieldSet & fsetFg,
                       const util::DateTime & validTimeOfXbFg,
                       std::vector<atlas::FieldSet> & fsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const std::vector<saber::SaberOuterBlockParametersWrapper> & params) {
  oops::Log::trace() << "SaberOuterBlockChain ctor starting" << std::endl;
  oops::Log::info() << "Info     : Creating outer blocks" << std::endl;

  // In addition to other configuration option pass model data information for vader
  // TODO(AS): check whether covarConf needs to be passed to the blocks (ideally not)
  oops::ModelData<MODEL> modelData{geom};
  eckit::LocalConfiguration modelDataConf;
  modelDataConf.set("model data", modelData.modelData());  // Key here is vader::configModelVarsKey
  eckit::LocalConfiguration outerBlockConf{covarConf};
  outerBlockConf.set("vader", modelDataConf);

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

  // Loop in reverse order
  for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
    boost::adaptors::reverse(params)) {
    // Initialize current outer variables and outer geometry data
    oops::Variables currentOuterVars = outerBlocks_.size() == 0 ?
                                       outerVars : outerBlocks_.back()->innerVars();
    const oops::GeometryData & currentOuterGeometryData = outerBlocks_.size() == 0 ?
                                       geom.generic() : outerBlocks_.back()->innerGeometryData();

    // Get outer block parameters
    const SaberBlockParametersBase & saberOuterBlockParams =
      saberOuterBlockParamWrapper.saberOuterBlockParameters;
    oops::Log::info() << "Info     : Creating outer block: "
                      << saberOuterBlockParams.saberBlockName.value() << std::endl;

    // Get active variables
    oops::Variables activeVars = getActiveVars(saberOuterBlockParams, currentOuterVars);

    // Create outer block
    outerBlocks_.emplace_back(SaberOuterBlockFactory::create(
                                                 currentOuterGeometryData,
                                                 currentOuterVars,
                                                 outerBlockConf,
                                                 saberOuterBlockParams,
                                                 fsetXb,
                                                 fsetFg,
                                                 validTimeOfXbFg));

    // Read and add model fields
    outerBlocks_.back()->read(geom, currentOuterVars, validTimeOfXbFg);

    if (saberOuterBlockParams.doCalibration()) {
      // Block calibration

      // Ensemble configuration
      eckit::LocalConfiguration ensembleConf
        = covarConf.getSubConfiguration("ensemble configuration");

      if (iterativeEnsembleLoading) {
        // Iterative calibration
         oops::Log::info() << "Info     : Iterative calibration" << std::endl;

        // Initialization
        outerBlocks_.back()->iterativeCalibrationInit();

        // Get ensemble size
        const size_t nens = ensembleConf.getInt("ensemble size");

        for (size_t ie = 0; ie < nens; ++ie) {
          // Read ensemble member
          atlas::FieldSet fset;
          readEnsembleMember(geom,
                             outerVars,
                             validTimeOfXbFg,
                             ensembleConf,
                             ie,
                             fset);

          // Apply outer blocks inverse (except last)
          this->leftInverseMultiplyExceptLast(fset);

          // Use FieldSet in the central block
          oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
          outerBlocks_.back()->iterativeCalibrationUpdate(fset);
        }
        // Finalization
        oops::Log::info() << "Info     : Finalization" << std::endl;
        outerBlocks_.back()->iterativeCalibrationFinal();
      } else {
        // Direct calibration
        oops::Log::info() << "Info     : Direct calibration" << std::endl;
        outerBlocks_.back()->directCalibration(fsetEns);
      }

      // Write calibration data
      oops::Log::info() << "Info     : Write calibration data" << std::endl;
      outerBlocks_.back()->write(geom, currentOuterVars, validTimeOfXbFg);
      outerBlocks_.back()->write();

      if (!iterativeEnsembleLoading) {
        // Left inverse multiplication on ensemble members
        oops::Log::info() << "Info     : Left inverse multiplication on ensemble members"
                          << std::endl;
        for (auto & fset : fsetEns) {
          if (outerBlocks_.back()->skipInverse()) {
            oops::Log::info()
                    << "Info     : Warning: left inverse multiplication skipped for block "
                    << outerBlocks_.back()->blockName() << std::endl;
          } else {
            outerBlocks_.back()->leftInverseMultiply(fset);
          }
        }
      }
    } else if (saberOuterBlockParams.doRead()) {
      // Read data
      oops::Log::info() << "Info     : Read data" << std::endl;
      outerBlocks_.back()->read();
    }

    // Inner geometry data and variables
    const oops::GeometryData & innerGeometryData =
      outerBlocks_.back()->innerGeometryData();
    const oops::Variables innerVars = outerBlocks_.back()->innerVars();

    // Check that active variables are present in either inner or outer variables, or both
    for (const auto & var : activeVars.variables()) {
      ASSERT(innerVars.has(var) || currentOuterVars.has(var));
    }

    // Get intersection of active variables and outer/inner variables
    oops::Variables activeOuterVars = currentOuterVars;
    activeOuterVars.intersection(activeVars);
    oops::Variables activeInnerVars = innerVars;
    activeInnerVars.intersection(activeVars);

    // Left inverse multiplication on xb and fg if inner and outer Geometry is different
    if (util::getGridUid(innerGeometryData.functionSpace())
      != util::getGridUid(currentOuterGeometryData.functionSpace())
      && saberOuterBlockParams.inverseVars.value().size() > 0) {
      oops::Log::info() << "Info     : Left inverse multiplication on xb and fg" << std::endl;
      // Share fields pointers
      atlas::FieldSet fsetXbInv;
      atlas::FieldSet fsetFgInv;
      for (const auto & var : saberOuterBlockParams.inverseVars.value().variables()) {
        fsetXbInv.add(fsetXb.field(var));
        fsetFgInv.add(fsetFg.field(var));
      }

      // Apply left inverse
      outerBlocks_.back()->leftInverseMultiply(fsetXbInv);
      outerBlocks_.back()->leftInverseMultiply(fsetFgInv);

      // Copy the fields back
      for (const auto & var : saberOuterBlockParams.inverseVars.value().variables()) {
        fsetXb.field(var) = fsetXbInv.field(var);
        fsetFg.field(var) = fsetFgInv.field(var);
      }
    }

    // Adjoint test
    if (covarConf.getBool("adjoint test")) {
      // Get tolerance
      const double localAdjointTolerance =
        saberOuterBlockParams.adjointTolerance.value().get_value_or(
        covarConf.getDouble("adjoint tolerance"));

      // Run test
      outerBlocks_.back()->adjointTest(geom.getComm(),
                                       currentOuterGeometryData,
                                       activeOuterVars,
                                       innerGeometryData,
                                       activeInnerVars,
                                       localAdjointTolerance);
    }
    // Inverse test
    const bool skipInverseTest = saberOuterBlockParams.skipInverseTest.value();
    if (covarConf.getBool("inverse test", false)) {
      oops::Log::info() << "Info     : Inverse test" << std::endl;
      if (skipInverseTest) {
        oops::Log::test() << "skipping inverse test for block "
                          << outerBlocks_.back()->blockName() << std::endl;
      } else {
        // Get inner and outer tolerances
        const double innerInverseTolerance = saberOuterBlockParams.innerInverseTolerance.value()
          .get_value_or(covarConf.getDouble("inverse tolerance"));
        const double outerInverseTolerance = saberOuterBlockParams.outerInverseTolerance.value()
          .get_value_or(covarConf.getDouble("inverse tolerance"));

        // Get inner and outer variables to compare
        oops::Variables innerVarsToCompare = saberOuterBlockParams.innerVariables.value()
          .get_value_or(activeInnerVars);
        oops::Variables outerVarsToCompare = saberOuterBlockParams.outerVariables.value()
          .get_value_or(activeOuterVars);

        // Run test
        outerBlocks_.back()->inverseTest(innerGeometryData,
                                         activeInnerVars,
                                         currentOuterGeometryData,
                                         activeOuterVars,
                                         innerVarsToCompare,
                                         outerVarsToCompare,
                                         innerInverseTolerance,
                                         outerInverseTolerance,
                                         geom.timeComm().rank());
      }
    }
  }
  oops::Log::trace() << "SaberOuterBlockChain ctor done" << std::endl;
}

}  // namespace saber
