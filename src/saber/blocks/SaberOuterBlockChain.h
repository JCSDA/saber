/*
 * (C) Copyright 2023- UCAR
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/interface/ModelData.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/DefaultCookbook.h"

#include "vader/vader.h"

namespace saber {

/// Chain of outer saber blocks (no central block). Can be used as the common
/// outer blocks for the hybrid covariance, outer blocks for static and ensemble
/// covariances, ensemble transform for the ensemble covariance.
class SaberOuterBlockChain {
 public:
  /// @brief Standard constructor using MODEL geometry
  template<typename MODEL>
  SaberOuterBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf,
                       const std::vector<SaberOuterBlockParametersWrapper> & params,
                       std::shared_ptr<oops::FieldSets> fsetEns = NULL,
                       const bool & centralDirectCalibration = false);
  /// @brief Simpler, limited constructor using only generic GeometryData
  SaberOuterBlockChain(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf,
                       const std::vector<SaberOuterBlockParametersWrapper> & params);

  ~SaberOuterBlockChain() = default;

  // Accessors
  // TODO(AS): this should be const (currently used to add ens transform blocks
  // to the outer blocks in SaberEnsembleBlockChain)
  std::vector<std::pair<std::shared_ptr<SaberOuterBlockBase>, bool>> & outerBlocks()
    {return outerBlocks_;}

  /// @brief Returns inner-most geometry data.
  const oops::GeometryData & innerGeometryData() const {
    if (outerBlocks_.back().second) {
      // Right-inverse mode
      return outerBlocks_.back().first->outerGeometryData();
    } else {
      // Direct mode
      return outerBlocks_.back().first->innerGeometryData();
    }
  }

  /// @brief Returns inner-most variables.
  const oops::Variables & innerVars() const {
    if (outerBlocks_.back().second) {
      // Right-inverse mode
      return outerBlocks_.back().first->outerVars();
    } else {
      // Direct modes
      return outerBlocks_.back().first->innerVars();
    }
  }

  /// @brief Forward multiplication by all outer blocks, 3D.
  void applyOuterBlocks(oops::FieldSet3D & fset3d) const {
    for (auto it = outerBlocks_.rbegin(); it != outerBlocks_.rend(); ++it) {
      if (it->second) {
        // Right-inverse mode
        it->first.get()->rightInverseMultiply(fset3d);
      } else {
        // Direct mode
        it->first.get()->multiply(fset3d);
      }
    }
  }

  /// @brief Adjoint multiplication by all outer blocks, 3D.
  void applyOuterBlocksAD(oops::FieldSet3D & fset3d) const {
    for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
      if (it->second) {
        // Right-inverse mode
        throw eckit::Exception("not implemented yet, but it should be", Here());
      } else {
        // Direct mode
        it->first.get()->multiplyAD(fset3d);
      }
    }
  }

  /// @brief Forward multiplication by all outer blocks, 4D.
  void applyOuterBlocks(oops::FieldSet4D & fset4d) const {
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      this->applyOuterBlocks(fset4d[jtime]);
    }
  }

  /// @brief Adjoint multiplication by all outer blocks, 4D.
  void applyOuterBlocksAD(oops::FieldSet4D & fset4d) const {
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      this->applyOuterBlocksAD(fset4d[jtime]);
    }
  }

  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the ones that haven't implemented inverse yet.
  void leftInverseMultiply(oops::FieldSet3D & fset) const {
    for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
      if (it->first.get()->skipInverse()) {
        oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                          << it->first.get()->blockName() << std::endl;
      } else {
        if (it->second) {
          // Right-inverse mode
          it->first.get()->multiply(fset);
        } else {
          // Direct mode
          it->first->leftInverseMultiply(fset);
        }
      }
    }
  }

  /// @brief Right inverse multiply (used in ensemble transform) by all outer blocks
  ///        except the ones that haven't implemented inverse yet.
  void rightInverseMultiply(oops::FieldSet3D & fset) const {
    for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
      if (it->first.get()->skipInverse()) {
        oops::Log::info() << "Warning: right inverse multiplication skipped for block "
                          << it->first.get()->blockName() << std::endl;
      } else {
        if (it->second) {
          // Right-inverse mode
          throw eckit::Exception("no right inverse available in right inverse mode", Here());
        } else {
          // Direct mode
          it->first.get()->rightInverseMultiply(fset);
        }
      }
    }
  }

 private:
  /// @brief Initialize outer block, and return tuple of current outer variables,
  ///        saber block parameters and active variables
  std::tuple<const SaberBlockParametersBase&,
             oops::Variables,
             oops::Variables>
     initBlock(const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper,
               const eckit::Configuration & outerBlockConf,
               const oops::GeometryData & outerGeometryData,
               const oops::Variables & outerVars,
               oops::FieldSet4D & fset4dXb,
               oops::FieldSet4D & fset4dFg);

  /// @brief Block calibration. Used in standard constructor.
  template<typename MODEL>
  void calibrateBlock(const eckit::Configuration & conf,
                      const oops::FieldSet4D & fset4dXb,
                      const oops::Geometry<MODEL> & geom,
                      const bool & validModelGeom,
                      const oops::Variables & outerVars,
                      const oops::Variables & currentOuterVars,
                      oops::FieldSets & fsetEns);

  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the last one and the ones that haven't implemented inverse yet.
  void leftInverseMultiplyExceptLast(oops::FieldSet3D & fset) const {
    // Outer blocks left inverse multiplication
    for (auto it = outerBlocks_.begin(); it != std::prev(outerBlocks_.end()); ++it) {
      if (it->first.get()->skipInverse()) {
        oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                          << it->first.get()->blockName() << std::endl;
      } else {
        if (it->second) {
          // Right-inverse mode
          it->first.get()->multiply(fset);
        } else {
          // Direct mode
          it->first.get()->leftInverseMultiply(fset);
        }
      }
    }
  }

  /// @brief Interpolate fields in background and first guess if inner and outer
  ///        geometryData are different. Used in constructors.
  void interpolateStates(
          const SaberBlockParametersBase & saberOuterBlockParams,
          const oops::GeometryData & outerGeometryData,
          oops::FieldSet4D & fset4dXb,
          oops::FieldSet4D & fset4dFg) const;

  /// @brief Inverse and adjoint test for last outer block. Used in constructors.
  void testLastOuterBlock(const eckit::Configuration & conf,
                          const SaberBlockParametersBase & saberOuterBlockParams,
                          const oops::GeometryData & outerGeometryData,
                          const oops::Variables & outerVars,
                          const oops::Variables & activeVars) const;

  /// @brief Vector of all outer blocks, paired with a boolean to indicate the application mode:
  /// - true: right inverse mode
  /// - false: direct mode
  /// TODO(AS): Need to expand this to create different outer blocks for different
  /// times for the 4D with multiple times on one MPI task.
  std::vector<std::pair<std::shared_ptr<SaberOuterBlockBase>, bool>> outerBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberOuterBlockChain::SaberOuterBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       const eckit::Configuration & conf,
                       const std::vector<saber::SaberOuterBlockParametersWrapper> & params,
                       std::shared_ptr<oops::FieldSets> fsetEns,
                       const bool & centralDirectCalibration) {
  oops::Log::trace() << "SaberOuterBlockChain ctor starting" << std::endl;
  oops::Log::info() << "Info     : Creating outer blocks" << std::endl;

  // In addition to other configuration option pass model data information for vader
  // TODO(AS): check whether conf needs to be passed to the blocks (ideally not)
  oops::ModelData<MODEL> modelData{geom};
  eckit::LocalConfiguration modelDataConf;
  modelDataConf.set("model data", modelData.modelData());  // Key here is vader::configModelVarsKey
  eckit::LocalConfiguration outerBlockConf{conf};
  outerBlockConf.set("vader", modelDataConf);

  // Copy vector of parameters
  std::vector<SaberOuterBlockParametersWrapper> innerParams = params;

  // Flag to check if the MODEL geometry is still valid
  bool validModelGeom = true;

  // Loop in reverse order
  for (int jb = params.size()-1; jb >= 0; --jb) {
    // Initialize current outer geometry data
    const oops::GeometryData & currentOuterGeometryData = outerBlocks_.size() == 0 ?
      geom.generic() : innerGeometryData();

    // Initialize outer block
    const auto[saberOuterBlockParams,
               currentOuterVars,
               activeVars]
            = initBlock(params[jb],
                          outerBlockConf,
                          currentOuterGeometryData,
                          outerVars,
                          fset4dXb,
                          fset4dFg);

    // Update MODEL geometry validity, by checking whether the inner geometry data returned by
    // the last outer block shares the same reference as its own outer geometry data
    validModelGeom = validModelGeom &&
      (&(innerGeometryData()) == &currentOuterGeometryData);

    // Read and add model fields
    outerBlocks_.back().first->read(geom, validModelGeom, currentOuterVars);

    // Remove element from inner parameters
    innerParams.pop_back();

    if (saberOuterBlockParams.doCalibration()) {
      // Block calibration
      calibrateBlock(conf,
                     fset4dXb,
                     geom,
                     validModelGeom,
                     outerVars,
                     currentOuterVars,
                     *fsetEns);
    } else if (saberOuterBlockParams.doRead()) {
      // Read data
      oops::Log::info() << "Info     : Read data" << std::endl;
      outerBlocks_.back().first->read();
    }

    if (saberOuterBlockParams.forceWrite.value()) {
      // Write data
      oops::Log::info() << "Info     : Write data" << std::endl;
      outerBlocks_.back().first->write(geom, validModelGeom, currentOuterVars);
      outerBlocks_.back().first->write();
    }

    if (!conf.getBool("iterative ensemble loading", false)) {
      // Check if the left inverse multiplication of this block on ensemble members if needed,
      // when either the central block or an inner outer block needs a direct calibration,
      // or if the final ensemble output is required
      bool applyLeftInverse = centralDirectCalibration;
      for (const auto & innerSaberOuterBlockParamWrapper : innerParams) {
        const SaberBlockParametersBase & innerSaberOuterBlockParams =
          innerSaberOuterBlockParamWrapper.saberOuterBlockParameters;
        applyLeftInverse = applyLeftInverse || innerSaberOuterBlockParams.doCalibration();
      }
      applyLeftInverse = applyLeftInverse || conf.has("output ensemble");

      if (applyLeftInverse) {
        // Left inverse multiplication on ensemble members
        oops::Log::info() << "Info     : Left inverse multiplication on ensemble members"
                        << std::endl;
        if (outerBlocks_.back().first->skipInverse()) {
            oops::Log::info()
                    << "Info     : Warning: left inverse multiplication skipped for block "
                    << outerBlocks_.back().first->blockName() << std::endl;
        } else if (fsetEns) {
          for (size_t jj = 0; jj < fsetEns->size(); ++jj) {
            outerBlocks_.back().first->leftInverseMultiply((*fsetEns)[jj]);
          }
        }
      }
    }

    // Left inverse multiplication on xb and fg if inner and outer Geometry is different
    interpolateStates(saberOuterBlockParams,
                      currentOuterGeometryData,
                      fset4dXb,
                      fset4dFg);

    // Adjoint test and inverse test
    testLastOuterBlock(conf,
                       saberOuterBlockParams,
                       currentOuterGeometryData,
                       currentOuterVars,
                       activeVars);
  }
  oops::Log::trace() << "SaberOuterBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberOuterBlockChain::calibrateBlock(
            const eckit::Configuration & conf,
            const oops::FieldSet4D & fset4dXb,
            const oops::Geometry<MODEL> & geom,
            const bool & validModelGeom,
            const oops::Variables & outerVars,
            const oops::Variables & currentOuterVars,
            oops::FieldSets & fsetEns) {
  oops::Log::trace() << "calibrateBlock starting" << std::endl;

  if (conf.getBool("iterative ensemble loading", false)) {
    // Iterative calibration
    oops::Log::info() << "Info     : Iterative calibration" << std::endl;

    // Initialization
    outerBlocks_.back().first->iterativeCalibrationInit();

    // Get ensemble size
    const size_t nens = getNensFromConfig(conf);

    for (size_t ie = 0; ie < nens; ++ie) {
      // Read ensemble member
      oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
      readEnsembleMember(geom,
                         outerVars,
                         conf,
                         ie,
                         fset);

      // Apply outer blocks inverse (except last)
      this->leftInverseMultiplyExceptLast(fset);

      // Use FieldSet in the central block
      oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
      outerBlocks_.back().first->iterativeCalibrationUpdate(fset);
    }

    // Finalization
    oops::Log::info() << "Info     : Finalization" << std::endl;
    outerBlocks_.back().first->iterativeCalibrationFinal();
  } else {
    // Direct calibration
    oops::Log::info() << "Info     : Direct calibration" << std::endl;
    outerBlocks_.back().first->directCalibration(fsetEns);
  }

  // Write calibration data
  oops::Log::info() << "Info     : Write calibration data" << std::endl;
  outerBlocks_.back().first->write(geom, validModelGeom, currentOuterVars);
  outerBlocks_.back().first->write();

  oops::Log::trace() << "calibrateBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
