/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/interface/ModelData.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

class ScaleParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ScaleParameters, Parameters)

 public:
  // Filter outer blockchain (optional)
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>> filterParams{
    "filter", this};

  // Use residual from filter
  oops::Parameter<bool> residualFromFilter{
    "use residual from filter", false, this};

  // Interpolator outer blockchain (optional)
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>> interpolatorParams{
    "interpolator", this};

  // Ensemble perturbations to read (optional)
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{
    "ensemble pert", this};

  // Ensemble perturbations to read on other geometry (optional)
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePertOtherGeom{
                        "ensemble pert on other geometry", this};

  // Output filtered perturbations (optional)
  oops::OptionalParameter<eckit::LocalConfiguration> output{
    "output", this};

  // Localization parametric blockchain (optional)
  oops::OptionalParameter<eckit::LocalConfiguration> localizationParams{
    "localization", this};

  // Localization at full resolution (even if an interpolator is present)
  oops::Parameter<bool> locAtFullRes{"localization at full resolution", false, this};
};


// -----------------------------------------------------------------------------

class ScaleData {
 public:
  explicit ScaleData(const ScaleParameters & params) :
    params_(params),
    internalInterpolation_(params_.interpolatorParams.value() && params_.locAtFullRes.value()),
    externalInterpolation_(params_.interpolatorParams.value() && !params_.locAtFullRes.value()) {}

  // Accessors
  const std::unique_ptr<SaberOuterBlockChain> & filter() const
    {return filter_;}
  std::unique_ptr<SaberOuterBlockChain> & filter()
    {return filter_;}
  const std::unique_ptr<SaberOuterBlockChain> & interpolator() const
    {return interpolator_;}
  std::unique_ptr<SaberOuterBlockChain> & interpolator()
    {return interpolator_;}
  const std::unique_ptr<SaberParametricBlockChain> & localization() const
    {return localization_;}
  std::unique_ptr<SaberParametricBlockChain> & localization()
    {return localization_;}
  const std::unique_ptr<oops::FieldSets> & ensemble() const
    {return ensemble_;}
  std::unique_ptr<oops::FieldSets> & ensemble()
    {return ensemble_;}
  const ScaleParameters & params() const
    {return params_;}
  const bool & internalInterpolation() const
    {return internalInterpolation_;}
  const bool & externalInterpolation() const
    {return externalInterpolation_;}

 private:
  /// @brief Filter outer block chain.
  std::unique_ptr<SaberOuterBlockChain> filter_;
  /// @brief Interpolator outer block chain.
  std::unique_ptr<SaberOuterBlockChain> interpolator_;
  /// @brief Localization parametric block chain.
  std::unique_ptr<SaberParametricBlockChain> localization_;
  /// @brief Ensemble used for this scale.
  std::unique_ptr<oops::FieldSets> ensemble_;
  /// @brief Scale parameters
  const ScaleParameters params_;
  /// @brief Internal interpolation
  bool internalInterpolation_;
  /// @brief External interpolation
  bool externalInterpolation_;
};

// -----------------------------------------------------------------------------

class InflationFieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InflationFieldParameters, Parameters)

 public:
  // ATLAS inflation file
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFileConf{"atlas file", this};
  // Model inflation file
  oops::OptionalParameter<eckit::LocalConfiguration> modelFileConf{"model file", this};
};

// -----------------------------------------------------------------------------

class SaberEnsembleBlockChainParameters: public ErrorCovarianceParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberEnsembleBlockChainParameters,
                           ErrorCovarianceParametersBase)

 public:
  // Outer blocks
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    saberOuterBlocksParams{"saber outer blocks", this};

  // Localization parameters
  oops::OptionalParameter<eckit::LocalConfiguration>
    localization{"localization", this};

  // Vector of scale-specific configurations
  oops::OptionalParameter<std::vector<ScaleParameters>> scales{"scales", this};

  // Recursive perturbations processing
  oops::Parameter<bool> recursivePertProcessing{"recursive perturbations processing", false, this};

  // Multi-scales strategy (separated or crossed)
  oops::OptionalParameter<std::string> strategy{"multiscale strategy", this};

  // Ensemble transform parameters
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    ensembleTransform{"ensemble transform", this};

  // Inflation fields
  oops::OptionalParameter<InflationFieldParameters> inflationField{"inflation field", this};

  // Inflation value
  oops::Parameter<double> inflationValue{"inflation value", 1.0, this};

  // Denominator for normalizing ensemble covariance (optional, default is to use ensemble size - 1)
  oops::OptionalParameter<double> denominatorForNormalizingEnsembleCovariance{
                        "denominator for normalizing ensemble covariance", this};

  // Ensemble
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleBase{"ensemble base", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePairs{"ensemble pairs", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePertOtherGeom{
                        "ensemble pert on other geometry", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensembleGeom{
                        "ensemble geometry", this};

  // Sub-ensembles size: it means that sets of members {0,...,subEnsSize-1},
  // {subEnsSize,...,2 subEnsSize-1}, etc. are distinct sub-ensembles. The mean of each sub-ensemble
  // is subtracted from the concerned members to compute perturbations.
  oops::OptionalParameter<size_t> subEnsSize{"sub-ensembles size", this};
};

/// Chain of outer (optional) and an ensemble "block".
class SaberEnsembleBlockChain : public SaberBlockChainBase {
 public:
  template<typename MODEL>
  SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                          const oops::Variables & outerVars,
                          oops::FieldSet4D & fset4dXb,
                          oops::FieldSet4D & fset4dFg,
                          const eckit::Configuration & conf);
  ~SaberEnsembleBlockChain() = default;

  /// @brief Randomize the increment according to this B matrix.
  void randomize(oops::FieldSet4D &) const;
  /// @brief Multiply the increment by this B matrix.
  void multiply(oops::FieldSet4D &) const;
  /// @brief Get this B matrix square-root control vector size.
  size_t ctlVecSize() const {return ctlVecSize_;}
  /// @brief Generate a random control vector.
  void randomCtlVec(atlas::Field &, const size_t &) const;
  /// @brief Multiply the control vector by this B matrix square-root.
  void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const;
  /// @brief Multiply the increment by this B matrix square-root adjoint.
  void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &) const;

  /// @brief Accessor to outer function space
  const atlas::FunctionSpace & outerFunctionSpace() const {return outerFunctionSpace_;}
  /// @brief Accessor to outer variables
  const oops::Variables & outerVariables() const {return outerVariables_;}

 private:
  /// @brief Space communicator
  const eckit::mpi::Comm & comm_;
  /// @brief Outer function space
  const atlas::FunctionSpace outerFunctionSpace_;
  /// @brief Outer variables
  const oops::Variables outerVariables_;
  /// @brief Outer blocks (optional).
  std::unique_ptr<SaberOuterBlockChain> outerBlockChain_;
  /// @brief Vector of data containers, one for each scale
  std::vector<ScaleData> scaleDataVec_;
  /// @brief Multiscales strategy
  std::string strategy_;
  /// @brief Control vector size.
  size_t ctlVecSize_;
  /// @brief Variables used in the ensemble covariance.
  /// TODO(AS): check whether this is needed or can be inferred from ensemble->
  oops::Variables vars_;
  int seed_ = 7;  // For reproducibility
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SaberEnsembleBlockChain::SaberEnsembleBlockChain(const oops::Geometry<MODEL> & geom,
                                                 const oops::Variables & outerVars,
                                                 oops::FieldSet4D & fset4dXb,
                                                 oops::FieldSet4D & fset4dFg,
                                                 const eckit::Configuration & conf)
  : comm_(geom.getComm()),
    outerFunctionSpace_(geom.functionSpace()),
    outerVariables_(outerVars),
    ctlVecSize_(0) {
  oops::Log::trace() << "SaberEnsembleBlockChain ctor starting" << std::endl;

  // Deserialize parameters and fill configuration with missing values
  SaberEnsembleBlockChainParameters params;
  params.deserialize(conf);
  eckit::LocalConfiguration fullConf;
  params.serialize(fullConf);

  // Extract ErrorCovarianceParametersBase from fullConf
  ErrorCovarianceParametersBase paramsBase;
  paramsBase.deserialize(fullConf);

  // Read ensemble (for non-iterative ensemble loading)
  std::unique_ptr<oops::FieldSets> ensemble = std::make_unique<oops::FieldSets>(
                                     readAndScaleEnsemble(
                                       geom,
                                       outerVars,
                                       fset4dXb.times(), fset4dXb.commTime(), fset4dXb.commEns(),
                                       fullConf));

  // Check that there is an ensemble of at least 2 members (or no member).
  if (ensemble->ens_size() == 1) {
    throw eckit::BadParameter("Ensemble for SaberEnsembleBlockChain has to have at least"
                              " two members (or no member)", Here());
  }

  // Remove specific means for sub-ensembles
  if (params.subEnsSize.value()) {
    // Sub-ensemble size
    const size_t subEnsSize = *params.subEnsSize.value();

    // Consistency checks
    ASSERT(subEnsSize > 1);
    ASSERT(ensemble->ens_size()%subEnsSize == 0);

    // Loop over sub-ensembles
    for (size_t jsub=0; jsub < ensemble->ens_size()/subEnsSize; ++jsub) {
      // Member index offset
      const size_t offset = jsub*subEnsSize;

      for (size_t it = 0; it < fset4dXb.size(); ++it) {
        // Initialize mean with first member of the sub-ensemble
        oops::FieldSet3D mean((*ensemble)(it, offset));

        // Accumule other members of the sub-ensemble
        for (size_t ie = 1; ie < subEnsSize; ++ie) {
          mean += (*ensemble)(it, offset+ie);
        }

        // Normalize mean
        const double rk = 1.0/static_cast<double>(subEnsSize);
        mean *= rk;

        // Subtract mean
        for (size_t ie = 0; ie < subEnsSize; ++ie) {
          (*ensemble)(it, offset+ie) -= mean;
        }
      }
    }
  }

  // Create outer blocks if needed
  if (params.saberOuterBlocksParams.value()) {
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(geom, outerVars,
                          fset4dXb, fset4dFg, fullConf,
                          *params.saberOuterBlocksParams.value());
  }

  // Outer variables and geometry for the ensemble covariance
  oops::Variables currentOuterVars = outerBlockChain_ ?
                                     outerBlockChain_->innerVars() : outerVars;
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                                   outerBlockChain_->innerGeometryData() : geom.generic();

  // Get active variables
  const oops::Variables activeVars = currentOuterVars;
  vars_ += activeVars;
  // Check that active variables are present in variables
  for (const auto & var : activeVars) {
    if (!currentOuterVars.has(var)) {
      throw eckit::UserError("Active variable " + var.name() + " is not present in "
                             "outer variables", Here());
    }
  }

  // Check consistency if ensemble was read on non-MODEL geometry
  if (fullConf.has("ensemble pert on other geometry")) {
    if (ensemble->ens_size() == 0) {
      throw eckit::BadParameter("No member loaded for perturbations on other geometry", Here());
    }
    const auto & currentFspace = currentOuterGeom.functionSpace();
    const auto & ensFspace = (*ensemble)[0].fieldSet()[0].functionspace();
    ASSERT(ensFspace.type() == currentFspace.type());
    ASSERT(util::getGridUid(ensFspace).compare(
               util::getGridUid(currentFspace)) == 0);
  }

  // Read inflation field
  const double inflationValue = params.inflationValue;
  oops::Log::info() << "Info     : Read inflation field" << std::endl;
  oops::FieldSet3D inflationField(fset4dXb[0].validTime(), geom.getComm());
  if (params.inflationField.value()) {
    const auto & inflationParams = *params.inflationField.value();
    // Read ATLAS inflation file
    if (inflationParams.atlasFileConf.value()) {
      eckit::LocalConfiguration inflationConf = *inflationParams.atlasFileConf.value();
      // Read file
      inflationField.read(currentOuterGeom.functionSpace(),
                          activeVars,
                          inflationConf);
      // Set name
      inflationField.name() = "inflation";

      // Print FieldSet norm
      oops::Log::test() << "Norm of input parameter inflation: "
                        << inflationField.norm(activeVars) << std::endl;
    }
    // Use model inflation file
    if (inflationParams.modelFileConf.value()) {
      eckit::LocalConfiguration inflationConf = *inflationParams.modelFileConf.value();
      // Copy file
      // Read fieldsets as increments
      // Create increment
      oops::Increment<MODEL> dx(geom, activeVars, fset4dXb[0].validTime());
      dx.read(inflationConf);
      oops::Log::test() << "Norm of input parameter inflation"
                        << ": " << dx.norm() << std::endl;
      inflationField.deepCopy(dx.fieldSet());
    }
  }

  if (ensemble->ens_size() > 0) {
    // Apply inflation on ensemble members
    oops::Log::info() << "Info     : Apply inflation on ensemble members" << std::endl;

    // Apply local inflation
    if (!inflationField.empty()) {
      *ensemble *= inflationField;
    }

    // Apply global inflation
    *ensemble *= inflationValue;

    // Ensemble transform
    if (params.ensembleTransform.value()) {
      const auto ensTransParams = *params.ensembleTransform.value();
      oops::Log::info() << "Info     : Found ensemble transform " << std::endl;
      std::unique_ptr<SaberOuterBlockChain> ensTransBlockChain =
             std::make_unique<SaberOuterBlockChain>(geom,
               currentOuterVars, fset4dXb, fset4dFg,
               paramsBase.toConfiguration(), ensTransParams);

      // Right inverse of ensemble transform on ensemble members
      oops::Log::info() << "Info     : Right inverse of ensemble transform on ensemble members"
                        << std::endl;
      for (size_t it = 0; it < ensemble->local_time_size(); ++it) {
        for (size_t ie = 0; ie < ensemble->local_ens_size(); ++ie) {
          ensTransBlockChain->rightInverseMultiply((*ensemble)(it, ie));
        }
      }

      // Add ensemble transform blocks to outer blocks
      // TODO(AS): refactor so there is no need for non-const accessor to outerBlocks
      // in SaberOuterBlockChain.
      oops::Log::info() << "Info     : Add ensemble transform blocks to outer blocks"
                        << std::endl;
      if (outerBlockChain_) {
        std::move(ensTransBlockChain->outerBlocks().begin(),
                  ensTransBlockChain->outerBlocks().end(),
                  std::back_inserter(outerBlockChain_->outerBlocks()));
        ensTransBlockChain->outerBlocks().erase(ensTransBlockChain->outerBlocks().begin(),
                                                ensTransBlockChain->outerBlocks().end());
      } else {
        outerBlockChain_ = std::move(ensTransBlockChain);
      }

      // Update outer variables
      currentOuterVars = outerBlockChain_->innerVars();

      // Check that the geometry is still the same
      ASSERT(&currentOuterGeom == &(outerBlockChain_->innerGeometryData()));
    }
  }

  // Get scales configurations for SDL
  std::vector<ScaleParameters> scalesParams;
  if (params.scales.value()) {
    // Get vector of configurations, one for each scale
    scalesParams = *params.scales.value();

    // Check localization presence
    for (const auto & scaleParams : scalesParams) {
      ASSERT(scaleParams.localizationParams.value());
    }

    // Get multi-scales strategy (separated or crossed)
    ASSERT(params.strategy.value());
    strategy_ = *params.strategy.value();
    ASSERT((strategy_ == "separated") || (strategy_ == "crossed"));
  } else {
    // No scale separation
    eckit::LocalConfiguration scaleConf;

    // Add localization if present
    if (params.localization.value()) {
      scaleConf.set("localization", *params.localization.value());
    }

    // Add configuration
    ScaleParameters scaleParams;
    scaleParams.deserialize(scaleConf);
    scalesParams.push_back(scaleParams);

    // Set strategy
    strategy_ = "separated";
  }

  // Loop over scales
  for (const auto & scaleParams : scalesParams) {
    oops::Log::info() << "Info     : Scale " << (scaleDataVec_.size()+1) << " setup" << std::endl;

    // Create data container for this scale
    ScaleData scaleData(scaleParams);

    // Consistency check 1: filter or not filter?
    if (scalesParams.size() > 1) {
      if (scalesParams[0].filterParams.value()) {
        // First scale include a filter: all scales should have one too, except the last one maybe
        if (scaleDataVec_.size() < scalesParams.size()-1) {
          ASSERT(scaleParams.filterParams.value());
        }
      } else {
        // First scale does not include a filter: all scales should read ensemble perturbations
        ASSERT(scaleParams.ensemblePert.value() || scaleParams.ensemblePertOtherGeom.value());
      }
    }

    // Consistency check 2: residual from filter is incompatible with interpolator
    if (scaleParams.residualFromFilter.value()) {
      ASSERT(!scaleParams.interpolatorParams.value());
    }

    // Get interpolator outer geometry data
    const oops::GeometryData & interpolatorOuterGeomData = outerBlockChain_ ?
      outerBlockChain_->innerGeometryData() : geom.generic();

    // Interpolator outer block chain
    if (scaleParams.interpolatorParams.value()) {
      oops::Log::info() << "Info     : Interpolator setup" << std::endl;

      // Initialize interpolator outer block chain
      scaleData.interpolator() = std::make_unique<SaberOuterBlockChain>(
        interpolatorOuterGeomData,
        currentOuterVars,
        fset4dXb,
        fset4dFg,
        paramsBase.toConfiguration(),
        *scaleParams.interpolatorParams.value());

      // Interpolator should not change the variables
      ASSERT(scaleData.interpolator()->innerVars() == currentOuterVars);
    }

    // Get filter outer geometry data
    const oops::GeometryData & filterOuterGeomData = scaleData.interpolator() ?
      scaleData.interpolator()->innerGeometryData() : interpolatorOuterGeomData;

    // Filter outer block chain
    if (scaleParams.filterParams.value()) {
      oops::Log::info() << "Info     : Filter setup" << std::endl;

      // Create configuration without tests for filters
      const ErrorCovarianceParametersBase defaultParamsBase;

      // Initialize filter outer block chain
      scaleData.filter() = std::make_unique<SaberOuterBlockChain>(
        filterOuterGeomData,
        currentOuterVars,
        fset4dXb,
        fset4dFg,
        defaultParamsBase.toConfiguration(),
        *scaleParams.filterParams.value());

      // Filter should not change the variables
      ASSERT(scaleData.filter()->innerVars() == currentOuterVars);

      // Check geometry data consistency
      ASSERT(util::getGridUid(scaleData.filter()->innerGeometryData().functionSpace())
        == util::getGridUid(interpolatorOuterGeomData.functionSpace()));
    }

    // Get localization geometry data
    const oops::GeometryData & localizationGeomData = scaleData.externalInterpolation() ?
       filterOuterGeomData : interpolatorOuterGeomData;

    // Localization
    if (scaleParams.localizationParams.value()) {
      oops::Log::info() << "Info     : Localization setup" << std::endl;

      // Merge localization configuration with full configuration (order of arguments matters!)
      const eckit::LocalConfiguration locMergedConf =
        util::mergeConfigs(*scaleParams.localizationParams.value(), paramsBase.toConfiguration());

      // The localization is a parametric block chain constructed with the same geometry
      // as the ensemble block chain by default. If the outer blocks or transform block
      // chain include a change of geometries, we need to use build the localization from
      // the current geometryData, using a generic constructor of the parametric block chain.

      // Check consistency of `geom` and the current geometry
      const auto & currentFspace = localizationGeomData.functionSpace();
      if (util::getGridUid(geom.functionSpace()) != util::getGridUid(currentFspace)) {
        oops::Log::info() << "Info     : Localization and ensemble are on different "
                             "functionSpaces, building localization with generic "
                             "constructor" << std::endl;
        // Note QUENCH could just build another geometry here and use the standard
        // constructor, but other models usually don't have this ability to create a
        // Geometry on any mesh.
        scaleData.localization() = std::make_unique<SaberParametricBlockChain>(localizationGeomData,
                                                                     currentOuterVars,
                                                                     fset4dXb,
                                                                     fset4dFg,
                                                                     locMergedConf);
      } else {
        oops::Log::info() << "Info     : Localization and ensemble are on same "
                             "functionSpaces, building localization with standard "
                             "constructor" << std::endl;
        scaleData.localization() = std::make_unique<SaberParametricBlockChain>(geom,
                                                                     currentOuterVars,
                                                                     fset4dXb,
                                                                     fset4dFg,
                                                                     locMergedConf);
      }
    } else {
      // No localization, only authorized if no scales are specified
      oops::Log::info() << "Info     : No localization" << std::endl;
      ASSERT(!params.scales.value());
    }

    // Add scale data
    scaleDataVec_.emplace_back(std::move(scaleData));
  }

  if (strategy_ == "crossed") {
    // Check that localization control vector size is the same for all scales
    const size_t ctlVecSize0 = scaleDataVec_[0].localization()->ctlVecSize();
    for (auto & scaleData : scaleDataVec_) {
      ASSERT(scaleData.localization()->ctlVecSize() == ctlVecSize0);
    }
  }

  // Prepare ensembles
  if (params.scales.value()) {
    if (scaleDataVec_[0].filter()) {
      // Split ensemble into scales
      oops::Log::info() << "Info     : Split ensemble into scales" << std::endl;

      // Create empty ensemble for each scale
      for (auto & scaleData : scaleDataVec_) {
        scaleData.ensemble() = std::make_unique<oops::FieldSets>(
                                 ensemble->times(),
                                 ensemble->commTime(),
                                 ensemble->members(),
                                 ensemble->commEns());
      }

      // Process members sequentially
      oops::Log::info() << "Info     : Ensemble member : ";
      for (size_t ie = 0; ie < ensemble->ens_size(); ++ie) {
        // Print ensemble member index
        oops::Log::info() << (ie+1);
        if (ie < ensemble->ens_size()-1) {
          oops::Log::info() << " ";
        } else {
          oops::Log::info() << std::endl;
        }

        // Loop over subwindows
        for (size_t it = 0; it < fset4dXb.size(); ++it) {
          // Initialize work perturbation xI from ensemble perturbation x0
          oops::FieldSet3D fsetI((*ensemble)(it, ie));

          // Initialize sum of filtered perturbations
          oops::FieldSet3D fsetSum(fsetI.validTime(), fsetI.commGeom());
          fsetSum.allocateOnly(fsetI.fieldSet());
          fsetSum.zero();
          oops::FieldSet4D fset4dDxSum(fsetSum);

          for (auto & scaleData : scaleDataVec_) {
            // Copy work perturbation x = xI
            oops::FieldSet3D fset(fsetI.validTime(), fsetI.commGeom());
            fset.deepCopy(fsetI.fieldSet());
            oops::FieldSet4D fset4dDx(fset);

            if (scaleData.filter()) {
              // Apply filter G on input x: x' = Gx
              scaleData.filter()->applyOuterBlocks(fset4dDx);

              if (scaleData.params().residualFromFilter.value()) {
                // Use filter complement: x' = (I-G)x
                fset4dDx[0] -= fsetI;
                fset4dDx[0] *= -1.0;

                if (params.recursivePertProcessing.value()) {
                  // Recursive processing: xI = xI - x'
                  fsetI -= fset4dDx[0];
                }

                // Increment sum with the latest x'
                fset4dDxSum += fset4dDx;
              } else {
                // If needed, interpolate copy of the filtered perturbation to full resolution
                std::unique_ptr<oops::FieldSet4D> fset4dDxUPtr{};
                oops::FieldSet4D * fset4dDxPtr = &fset4dDx;
                if (scaleData.interpolator()) {
                  fset4dDxUPtr = std::unique_ptr<oops::FieldSet4D>(new oops::FieldSet4D(fset4dDx));
                  fset4dDxPtr = fset4dDxUPtr.get();
                  scaleData.interpolator()->applyOuterBlocks(*fset4dDxPtr);
                }

                if (params.recursivePertProcessing.value()) {
                  // Recursive processing: xI = xI - x'
                  fsetI -= (*fset4dDxPtr)[0];
                }

                // Increment sum with the latest x'
                fset4dDxSum += *fset4dDxPtr;
              }
            } else {
              // No filter on the last scale, use residual increment: x' = x0 - sum{previous x'}
              fset4dDx[0].zero();
              fset4dDx[0] += (*ensemble)(it, ie);
              fset4dDx[0] -= fset4dDxSum[0];
            }

            // Copy perturbation into ensemble
            scaleData.ensemble()->emplace_back(it, ie, fset4dDx[0]);

            if (&scaleData == &scaleDataVec_.back()) {
              // Remove member from initial ensemble as there are not needed anymore
              ensemble->clear(it, ie);
            }
          }
        }
      }
    } else {
      // Read ensemble perturbations
      for (auto & scaleData : scaleDataVec_) {
        if (scaleData.params().ensemblePert.value()) {
          // Check geometry consistency
          if (outerBlockChain_) {
            ASSERT(util::getGridUid(outerBlockChain_->innerGeometryData().functionSpace())
              == util::getGridUid(geom.functionSpace()));
          }

          // No interpolator allowed
          ASSERT(!scaleData.interpolator());

          // Read ensemble using model reader
          scaleData.ensemble() = std::make_unique<oops::FieldSets>(readEnsemble(
                                   geom,
                                   scaleData.localization()->outerVariables(),
                                   fset4dXb.times(), fset4dXb.commTime(), fset4dXb.commEns(),
                                   scaleData.params().toConfiguration()));
        } else if (scaleData.params().ensemblePertOtherGeom.value()) {
          // Get current geometry
          const oops::GeometryData & ensGeom = scaleData.interpolator() ?
            scaleData.interpolator()->innerGeometryData() : outerBlockChain_ ?
            outerBlockChain_->innerGeometryData() : geom.generic();

          // Read ensemble
          scaleData.ensemble() = std::make_unique<oops::FieldSets>(ensGeom.functionSpace(),
                                  currentOuterVars, fset4dXb.times(),
                                  *scaleData.params().ensemblePertOtherGeom.value(),
                                  ensGeom.comm(), fset4dXb.commTime());
        }
      }
    }
  } else {
    // No scales, move ensemble pointer
    scaleDataVec_[0].ensemble() = std::move(ensemble);
  }

  // Perturbations output
  for (const auto & scaleData : scaleDataVec_) {
    if (scaleData.params().output.value()) {
      // Get number of filtered perturbations to write
      size_t ne = scaleData.ensemble()->ens_size();
      if (scaleData.params().output.value()->getBool("only first perturbation", false)) {
        ne = 1;
      }

      // Write filtered perturbations
      const eckit::LocalConfiguration oConf = *scaleData.params().output.value();
      if (oConf.has("generic write")) {
        const eckit::LocalConfiguration gConf = oConf.getSubConfiguration("generic write");
        for (size_t ie = 0; ie < ne; ++ie) {
          eckit::LocalConfiguration gConfMem(gConf);
          util::setMember(gConfMem, ie+1);
          util::writeFieldSet(comm_, gConfMem, (*scaleData.ensemble())(0, ie).fieldSet());
        }
      }
      if (oConf.has("model write")) {
        const eckit::LocalConfiguration mConf = oConf.getSubConfiguration("model write");
        for (size_t ie = 0; ie < ne; ++ie) {
          // Should be on the model geometry!
          auto pert = oops::Increment<MODEL>(geom,
                                             scaleData.localization()->outerVariables(),
                                             fset4dXb.times()[0]);
          pert.zero();
          pert.fromFieldSet((*scaleData.ensemble())(0, ie).fieldSet());

          eckit::LocalConfiguration mConfMem(mConf);
          util::setMember(mConfMem, ie+1);
          pert.write(mConfMem);
        }
      }
    }
  }

  // Get control vector size
  if (scaleDataVec_[0].localization()) {
    // Check that all scales have a localization
    for (const auto & scaleData : scaleDataVec_) {
      ASSERT(scaleData.localization());
    }

    // Compute control vector size
    if (strategy_ == "separated") {
      // Separated strategy
      for (const auto & scaleData : scaleDataVec_) {
        ctlVecSize_ += scaleData.ensemble()->ens_size()*scaleData.localization()->ctlVecSize();
      }
    } else if (strategy_ == "crossed") {
      // Crossed strategy
      ctlVecSize_ = scaleDataVec_[0].ensemble()->ens_size()
        *scaleDataVec_[0].localization()->ctlVecSize();

      // Check that all the scales have the same control vector size
      for (const auto & scaleData : scaleDataVec_) {
        ASSERT(scaleData.localization()->ctlVecSize() ==
          scaleDataVec_[0].localization()->ctlVecSize());
      }
    }
  } else {
    // Without localization
    // Only one scale allowed
    ASSERT(scaleDataVec_.size() == 1);

    // Control vector size = number of members
    ctlVecSize_ = scaleDataVec_[0].ensemble()->ens_size();
  }

  // Adjoint test
  // TODO(AS): this is now a copy of the test in SaberCentralBlock; needs to be generalized.
  // (Perhaps the adjoint[s] test can be moved to SaberBlockChainBase.
  if (fullConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance = params.adjointTolerance.value();

    // Create random FieldSets
    oops::FieldSet4D fset4d1(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d1.size(); ++jtime) {
      fset4d1[jtime].randomInit(currentOuterGeom.functionSpace(), activeVars);
    }
    oops::FieldSet4D fset4d2(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d2.size(); ++jtime) {
      fset4d2[jtime].randomInit(currentOuterGeom.functionSpace(), activeVars);
    }
    // Copy FieldSets
    const oops::FieldSet4D fset4d1Save = oops::copyFieldSet4D(fset4d1);
    const oops::FieldSet4D fset4d2Save = oops::copyFieldSet4D(fset4d2);

    // Apply forward multiplication only (self-adjointness test)
    // TODO(AS): need to change this to only call it for the ensemble part, not outer blocks!
    this->multiply(fset4d1);
    this->multiply(fset4d2);

    // Compute adjoint test
    const double dp1 = fset4d1.dot_product_with(fset4d2Save, activeVars);
    const double dp2 = fset4d2.dot_product_with(fset4d1Save, activeVars);
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : adjoint tolerance = "
                      << localAdjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for block Ensemble";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localAdjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Adjoint test failure for block Ensemble", Here());
    }
  }

  // Square-root test
  if (fullConf.getBool("square-root test")) {
    // Get tolerance
    const double localSqrtTolerance = params.sqrtTolerance.value();

    // Create FieldSet
    oops::FieldSet4D fset4d(fset4dXb.times(), fset4dXb.commTime(), currentOuterGeom.comm());
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      fset4d[jtime].randomInit(outerFunctionSpace_, outerVariables_);
    }

    // Copy FieldSet
    oops::FieldSet4D fset4dSave = oops::copyFieldSet4D(fset4d);

    // Create control vector
    oops::Log::info() << "Info     : Control vector size for block Ensemble: "
                      << ctlVecSize_ << std::endl;
    atlas::Field cv = atlas::Field("genericCtlVec",
                                   atlas::array::make_datatype<double>(),
                                   atlas::array::make_shape(ctlVecSize_));
    util::NormalDistribution<double> dist(ctlVecSize_, 0.0, 1.0, seed_);
    std::vector<double> randVec;
    for (size_t jnode = 0; jnode < ctlVecSize_; ++jnode) {
      randVec.push_back(dist[jnode]);
    }
    if (!scaleDataVec_[0].localization()) {
      currentOuterGeom.comm().broadcast(randVec, 0);
    }
    auto view = atlas::array::make_view<double, 1>(cv);
    for (size_t jnode = 0; jnode < ctlVecSize_; ++jnode) {
      view(jnode) = randVec[jnode];
    }

    // Copy control vector
    atlas::Field ctlVecSave = atlas::Field("genericCtlVec",
                                           atlas::array::make_datatype<double>(),
                                           atlas::array::make_shape(ctlVecSize_));
    auto viewSave = atlas::array::make_view<double, 1>(ctlVecSave);
    viewSave.assign(view);

    // Apply square-root multiplication
    this->multiplySqrt(ctlVecSave, fset4d, 0);

    // Apply square-root adjoint multiplication
    this->multiplySqrtAD(fset4dSave, cv, 0);

    // Compute adjoint test
    const double dp1 = fset4d.dot_product_with(fset4dSave, activeVars);
    double dp2 = 0.0;
    for (size_t jnode = 0; jnode < ctlVecSize_; ++jnode) {
      dp2 += view(jnode)*viewSave(jnode);
    }
    if (scaleDataVec_[0].localization()) {
      currentOuterGeom.comm().allReduceInPlace(dp2, eckit::mpi::sum());
      fset4d.commTime().allReduceInPlace(dp2, eckit::mpi::sum());
    }
    oops::Log::info() << std::setprecision(16) << "Info     : Square-root test: y^t (Ax) = " << dp1
                      << ": x^t (A^t y) = " << dp2 << " : square-root tolerance = "
                      << localSqrtTolerance << std::endl;
    const bool adjComparison = (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < localSqrtTolerance);

    // Apply square-root multiplication
    this->multiplySqrt(cv, fset4d, 0);

    // Apply full multiplication
    this->multiply(fset4dSave);

    // Check that the fieldsets are similar within tolerance
    bool sqrtComparison = true;
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      sqrtComparison = sqrtComparison && fset4d[jtime].compare_with(fset4dSave[jtime],
        localSqrtTolerance, util::ToleranceType::relative);
    }
    if (sqrtComparison) {
      oops::Log::info() << "Info     : Square-root test passed: U U^t x == B x" << std::endl;
    } else {
      oops::Log::info() << "Info     : Square-root test failed: U U^t x != B x" << std::endl;
    }

    // Print results
    oops::Log::test() << "Square-root test for block Ensemble";
    if (adjComparison && sqrtComparison) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Square-root test failure for block Ensemble", Here());
    }
  }

  oops::Log::trace() << "SaberEnsembleBlockChain ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
