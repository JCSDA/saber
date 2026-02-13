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
#include "saber/blocks/SaberCentralBlock.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

class SaberParametricBlockChainParameters: public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberParametricBlockChainParameters,
                           ErrorCovarianceParametersBase)

 public:
  // Central and outer blocks
  oops::RequiredParameter<SaberCentralBlockParameters>
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
  oops::Variables initCentralBlock(const oops::GeometryData & outerGeom,
                                   const eckit::Configuration & conf,
                                   const SaberCentralBlockParameters & saberCentralBlockParams,
                                   const oops::FieldSet4D & fset4dXb,
                                   const oops::FieldSet4D & fset4dFg);

  /// @brief Run adjoint and square-root tests on central block. Used in constructors.
  void testCentralBlock(const eckit::Configuration & conf) const;

  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  std::shared_ptr<SaberOuterBlockChain> outerBlockChain_;
  bool crossTimeCov_;
  std::unique_ptr<SaberCentralBlock> centralBlock_;
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
  SaberCentralBlockParameters saberCentralBlockParams = params.saberCentralBlockParams;
  const bool centralDirectCalibration = saberCentralBlockParams.doCalibration();

  // Read ensemble (for non-iterative ensemble loading)
  std::shared_ptr<oops::FieldSets> fsetEns = std::make_shared<oops::FieldSets>(readEnsemble(geom,
                                         outerVars,
                                         fset4dXb.times(), fset4dXb.commTime(), fset4dXb.commEns(),
                                         fullConf));

  // If needed create outer block chain
  if (params.saberOuterBlocksParams.value()) {
    outerBlockChain_ = std::make_shared<SaberOuterBlockChain>(geom, outerVariables_,
                          fset4dXb, fset4dFg, fullConf,
                          *params.saberOuterBlocksParams.value(),
                          fsetEns, centralDirectCalibration);
  }

  // Set outer geometry data for central block
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : geom.generic();

  // Create central block
  oops::Log::info() << "Info     : Creating central block: " << std::endl;

  const auto currentOuterVars = initCentralBlock(currentOuterGeom,
                                                 fullConf,
                                                 saberCentralBlockParams,
                                                 fset4dXb,
                                                 fset4dFg);

  // Read and add model fields
  centralBlock_->read(geom, currentOuterVars);

  if (centralBlock_->doCalibration()) {
    // Calibration
    centralBlock_->calibrateBlock(geom,
                                  outerVariables_,
                                  fset4dXb,
                                  fset4dFg,
                                  fullConf,
                                  outerBlockChain_,
                                  fsetEns);
  }

  if (centralBlock_->doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    centralBlock_->read();
  }

  if (centralBlock_->forceWrite() || centralBlock_->doCalibration()) {
    // Write data
    oops::Log::info() << "Info     : Write data" << std::endl;
    centralBlock_->write(geom);
    centralBlock_->write();
  }

  // Test central block
  testCentralBlock(fullConf);

  oops::Log::trace() << "SaberParametricBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
