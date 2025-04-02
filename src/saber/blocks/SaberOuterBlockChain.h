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
                       oops::FieldSets & fsetEns,
                       const eckit::LocalConfiguration & covarConf,
                       const std::vector<saber::SaberOuterBlockParametersWrapper> & params);
  /// @brief Simpler, limited constructor using only generic GeometryData
  SaberOuterBlockChain(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
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
  void applyOuterBlocks(oops::FieldSet4D & fset4d) const {
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      for (auto it = outerBlocks_.rbegin(); it != outerBlocks_.rend(); ++it) {
        it->get()->multiply(fset4d[jtime]);
      }
    }
  }

  /// @brief Adjoint multiplication by all outer blocks.
  void applyOuterBlocksAD(oops::FieldSet4D & fset4d) const {
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
        it->get()->multiplyAD(fset4d[jtime]);
      }
    }
  }

  /// @brief Adjoint multiplication or filter to outer blocks.
  void applyOuterBlocksFilter(oops::FieldSet4D & fset) const {
    for (size_t jtime = 0; jtime < fset.size(); ++jtime) {
      for (const auto & outerBlocks : outerBlocks_) {
        if (outerBlocks->filterMode()) {
          outerBlocks->leftInverseMultiply(fset[jtime]);
        } else {
          outerBlocks->multiplyAD(fset[jtime]);
        }
      }
    }
  }

  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the ones that haven't implemented inverse yet.
  void leftInverseMultiply(oops::FieldSet3D & fset) const {
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
  /// @brief Initialize outer block, and return tuple of current outer variables,
  ///        saber block parameters and active variables
  std::tuple<const SaberBlockParametersBase&,
             oops::Variables,
             oops::Variables>
     initBlock(const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper,
               const eckit::LocalConfiguration & outerBlockConf,
               const oops::GeometryData & outerGeometryData,
               const oops::Variables & outerVars,
               oops::FieldSet4D & fset4dXb,
               oops::FieldSet4D & fset4dFg);

  /// @brief Block calibration. Used in standard constructor.
  template<typename MODEL>
  void calibrateBlock(const eckit::LocalConfiguration & covarConf,
                      const oops::FieldSet4D & fset4dXb,
                      const oops::Geometry<MODEL> & geom,
                      const oops::Variables & outerVars,
                      oops::FieldSets & fsetEns);

  /// @brief Left inverse multiply (used in calibration) by all outer blocks
  ///        except the last one and the ones that haven't implemented inverse yet.
  void leftInverseMultiplyExceptLast(oops::FieldSet3D & fset) const {
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

  /// @brief Get inner geometry data and variables, and check consistency with
  ///        active variables. Used in constructors.
  std::tuple<const oops::GeometryData &, const oops::Variables &>
        getInnerObjects(const oops::Variables & activeVars,
                        const oops::Variables & outerVars) const;

  /// @brief Interpolate fields in background and first guess if inner and outer
  ///        geometryData are different. Used in constructors.
  void interpolateStates(
          const SaberBlockParametersBase & saberOuterBlockParams,
          const oops::GeometryData & outerGeometryData,
          const oops::GeometryData & innerGeometryData,
          oops::FieldSet4D & fset4dXb,
          oops::FieldSet4D & fset4dFg) const;

  /// @brief Inverse and adjoint test for last outer block. Used in constructors.
  void testLastOuterBlock(const eckit::LocalConfiguration & covarConf,
                          const SaberBlockParametersBase & saberOuterBlockParams,
                          const oops::GeometryData & outerGeometryData,
                          const oops::Variables & outerVars,
                          const oops::GeometryData & innerGeometryData,
                          const oops::Variables & innerVars,
                          const oops::Variables & activeVars) const;

  /// @brief Vector of all outer blocks.
  /// TODO(AS): Need to expand this to create different outer blocks for different
  /// times for the 4D with multiple times on one MPI task.
  std::vector<std::unique_ptr<SaberOuterBlockBase>> outerBlocks_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberOuterBlockChain::SaberOuterBlockChain(const oops::Geometry<MODEL> & geom,
                       const oops::Variables & outerVars,
                       oops::FieldSet4D & fset4dXb,
                       oops::FieldSet4D & fset4dFg,
                       oops::FieldSets & fsetEns,
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

  // Loop in reverse order
  for (const SaberOuterBlockParametersWrapper & saberOuterBlockParamWrapper :
    boost::adaptors::reverse(params)) {
    // Initialize current outer geometry data
    const oops::GeometryData & currentOuterGeometryData = outerBlocks_.size() == 0 ?
                                       geom.generic() : outerBlocks_.back()->innerGeometryData();

    // Initialize outer block
    const auto[saberOuterBlockParams,
               currentOuterVars,
               activeVars]
            = initBlock(saberOuterBlockParamWrapper,
                          outerBlockConf,
                          currentOuterGeometryData,
                          outerVars,
                          fset4dXb,
                          fset4dFg);

    // Read and add model fields
    outerBlocks_.back()->read(geom, currentOuterVars);

    if (saberOuterBlockParams.doCalibration()) {
      // Block calibration
      calibrateBlock(covarConf,
                     fset4dXb,
                     geom,
                     currentOuterVars,
                     fsetEns);
    } else if (saberOuterBlockParams.doRead()) {
      // Read data
      oops::Log::info() << "Info     : Read data" << std::endl;
      outerBlocks_.back()->read();

      if (saberOuterBlockParams.forceWrite.value()) {
        // Write data
        oops::Log::info() << "Info     : Write data" << std::endl;
        outerBlocks_.back()->write(geom, outerVars);
        outerBlocks_.back()->write();
      }
    }

    // Inner geometry data and variables & consistency check with active variables
    auto[innerGeometryData, innerVars] = getInnerObjects(activeVars, currentOuterVars);

    // Left inverse multiplication on xb and fg if inner and outer Geometry is different
    interpolateStates(saberOuterBlockParams,
                      currentOuterGeometryData,
                      innerGeometryData,
                      fset4dXb,
                      fset4dFg);

    // Adjoint test and inverse test
    testLastOuterBlock(covarConf,
                       saberOuterBlockParams,
                       currentOuterGeometryData,
                       currentOuterVars,
                       innerGeometryData,
                       innerVars,
                       activeVars);
  }
  oops::Log::trace() << "SaberOuterBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberOuterBlockChain::calibrateBlock(
            const eckit::LocalConfiguration & covarConf,
            const oops::FieldSet4D & fset4dXb,
            const oops::Geometry<MODEL> & geom,
            const oops::Variables & outerVars,
            oops::FieldSets & fsetEns) {
  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = covarConf.getBool("iterative ensemble loading");

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
      oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
      readEnsembleMember(geom,
                         outerVars,
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
  outerBlocks_.back()->write(geom, outerVars);
  outerBlocks_.back()->write();

  if (!iterativeEnsembleLoading) {
    // Left inverse multiplication on ensemble members
    oops::Log::info() << "Info     : Left inverse multiplication on ensemble members"
                      << std::endl;
    if (outerBlocks_.back()->skipInverse()) {
        oops::Log::info()
                << "Info     : Warning: left inverse multiplication skipped for block "
                << outerBlocks_.back()->blockName() << std::endl;
    } else {
      for (size_t jj = 0; jj < fsetEns.size(); ++jj) {
        outerBlocks_.back()->leftInverseMultiply(fsetEns[jj]);
      }
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace saber
