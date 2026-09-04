/*
 * (C) Crown Copyright 2024 Met Office
 * (C) Copyright 2025 Meteorologisk Institutt
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/blocks/SaberCentralBlock.h"

#include <algorithm>

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Random.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {

// -----------------------------------------------------------------------------

bool SaberCentralBlockParameters::doCalibration() const {
  oops::Log::trace() << "SaberCentralBlockParameters::doCalibration starting" << std::endl;

  if (this->singleBlock.value()) {
    // Single block
    return this->singleBlock.value()->doCalibration();
  }
  if (this->groups.value()) {
    // Multiple blocks
    for (const auto & groupParams : this->groups.value().get()) {
      if (groupParams.centralBlockParams().doCalibration()) {
        return true;
      }
    }
  }
  return false;
}

// -----------------------------------------------------------------------------

bool SaberCentralBlockParameters::doRead() const {
  oops::Log::trace() << "SaberCentralBlockParameters::doRead starting" << std::endl;

  if (this->singleBlock.value()) {
    // Single block
    return this->singleBlock.value()->doRead();
  }
  if (this->groups.value()) {
    // Multiple blocks
    for (const auto & groupParams : this->groups.value().get()) {
      if (groupParams.centralBlockParams().doRead()) {
        return true;
      }
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
// TODO(anyone): needs cleanup

oops::Variables SaberCentralBlockParameters::getActiveVars(
  const oops::Variables & defaultVars) const {
  // Get groups configurations

  if (this->singleBlock.value()) {
    // Single block
    SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
    saberCentralBlockParamsWrapper.deserialize(this->toConfiguration());
    return saberCentralBlockParamsWrapper.blockParams().getActiveVars(defaultVars);
  } else {
    // Multiple groups
    oops::Variables activeVars_nomd;
    for (const auto & groupParams : *this->groups.value()) {
      const auto & params = groupParams.centralBlockParams();
      if (params.mandatoryActiveVars().size() == 0) {
        // No mandatory active variables for this block
        activeVars_nomd += params.activeVars.value().get_value_or(defaultVars);
      } else {
        // Block with mandatory active variables
        activeVars_nomd += params.activeVars.value().get_value_or(params.mandatoryActiveVars());
        ASSERT(params.mandatoryActiveVars() <= activeVars_nomd);
      }
    }

    // Copy the variables that exist in defaultVars from defaultVars (they have metadata
    // associated with them)
    oops::Variables activeVars;
    for (auto & var : activeVars_nomd) {
      if (defaultVars.has(var.name())) {
        activeVars.push_back(defaultVars[var.name()]);
      } else {
        activeVars.push_back(var);
      }
    }

    return activeVars;
  }
}

// -----------------------------------------------------------------------------

SaberCentralBlock::SaberCentralBlock(const oops::GeometryData & outerGeom,
                                     const oops::Variables & outerVars,
                                     const eckit::Configuration & covarConf,
                                     const SaberCentralBlockParameters & params,
                                     const oops::FieldSet3D & xb,
                                     const oops::FieldSet3D & fg)
  : geometryData_(outerGeom), outerVars_(outerVars), validTime_(xb.validTime()), params_(params),
    strategy_(params.strategy) {
  oops::Log::trace() << "SaberCentralBlock constructor starting" << std::endl;

  // Check multivariate strategy:
  // - Univariate: localization of each group is applied to each variable of the group.
  // - Duplicated: the localization of each group is to the sum of all the fields of the group,
  //   the result is split into the different fields.
  // - Duplicated and weighted: the localization of each group is applied to each variable,
  //   but variables of the same group are combined with user-specified weights.
  // - Crossed: the localization of each group is to the sum of all the fields of the group,
  //   the result is split into the different fields. All the groups share the same control
  //   vector: square-root formulation is necessary.
  ASSERT(strategy_ == "single" ||
         strategy_ == "univariate" ||
         strategy_ == "duplicated" ||
         strategy_ == "duplicated and weighted" ||
         strategy_ == "crossed");

  oops::Log::info() << "Info     : SaberCentralBlock using multivariate strategy: "
                    << strategy_ << std::endl;
  // Check if it's a single group
  if (params.singleBlock.value() && (strategy_ != "single")) {
    throw eckit::UserError("SaberCentralBlock: single block can only be used with the "
                           "'single' strategy.", Here());
  }
  if (strategy_ == "single" && params.groups.value() &&
             (params.groups.value().get().size() > 1)) {
    throw eckit::UserError("SaberCentralBlock: 'single' strategy can only be used with "
                           "a single block.", Here());
  }
  if (params.singleBlock.value() && params.groups.value()) {
    throw eckit::UserError("SaberCentralBlock: please specify either a single block or multiple "
                           "groups, not both.", Here());
  }

  if (params.singleBlock.value()) {
    // Single block case
    oops::Log::info() << "Info     : Creating single central block." << std::endl;

    // Number of groups
    ngroup_ = 1;

    // Append group parameters
    groupInputVars_.push_back(outerVars_);
    groupNames_.push_back("single block");
    groupInnerVars_.push_back(outerVars_);
    doCalibration_.push_back(params.singleBlock.value()->doCalibration());
    doRead_.push_back(params.singleBlock.value()->doRead());
    forceWrite_.push_back(params.singleBlock.value()->forceWrite.value());

    // Append empty outer block chain
    std::unique_ptr<SaberOuterBlockChain> emptyOuterBlockChainPtr;
    groupOuterBlockChains_.push_back(std::move(emptyOuterBlockChainPtr));

    // Append group central block
    groupCentralBlocks_.push_back(SaberCentralBlockFactory::create(outerGeom,
                                                       groupInnerVars_.back(),
                                                       covarConf,
                                                       *(params.singleBlock.value()),
                                                       xb,
                                                       fg));
  } else {
    // Number of groups
    ngroup_ = params.groups.value()->size();

    // Loop over groups
    for (const auto & groupParams : *params.groups.value()) {
      // Get group variables names
      oops::Log::info() << "Info     : Creating multivariate central block group: "
                        << groupParams.groupName.value() << std::endl;
      const oops::Variables varNames = groupParams.variables.value();

      // Define group variables
      oops::Variables groupVars;
      for (const auto & varName : varNames) {
        groupVars.push_back(outerVars_[varName.name()]);
      }

      // Check group variables consistency
      size_t levelsCheck = groupVars[0].getLevels();
      oops::Variable refVar = groupVars[0];
      for (const auto & var : groupVars) {
        // Get number of levels
        const size_t varLevels = var.getLevels();

        if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
          // All the fields of the group should have the same number of levels.
          ASSERT(levelsCheck == varLevels);
        } else if ((strategy_ == "duplicated") || (strategy_ == "crossed")) {
          // All the fields of the group should have the same number of levels,
          // or only one level if 3D and 2D fields are mixed (in this case, the reference
          // field should be 3D).
          if (varLevels > 1) {
            if (levelsCheck == 1) {
              levelsCheck = varLevels;
              refVar = var;
            } else {
              ASSERT(levelsCheck == varLevels);
            }
          }
        }
      }

      // Get variables on extra level
      const std::vector<std::string> & varsOnExtraLevels = groupParams.varsOnExtraLevels.value();
      ASSERT((varsOnExtraLevels.size() == 0) || (strategy_ == "duplicated"));

      // Initialize flags
      bool addFirstLevel = false;
      bool addLastLevel = false;

      for (auto & var : groupVars) {
        // First level
        int firstLevel = 0;

        // Check wheter an extra level is needed for this variable
        const bool extraLevel = (std::find(varsOnExtraLevels.begin(), varsOnExtraLevels.end(),
        var.name()) != varsOnExtraLevels.end());

        if ((refVar.getLevels() > 1) && (var.getLevels() == 1)) {
          // This is a 2D variable in a group containing 3D variables

          // Get field
          const auto field = xb[var.name()];

          // Find first level
          const std::string nearest3dLevel = field.metadata().getString("nearest 3d level");
          ASSERT((nearest3dLevel == "top") || (nearest3dLevel == "bottom"));
          if ((outerGeom.levelsAreTopDown() && (nearest3dLevel == "top")) ||
            (!outerGeom.levelsAreTopDown() && (nearest3dLevel == "bottom"))) {
            firstLevel = extraLevel ? -1 : 0;
          } else {
            firstLevel = extraLevel ? refVar.getLevels() : refVar.getLevels()-1;
          }

          // Update flags to add levels
          if (firstLevel == -1) {
            addFirstLevel = true;
          }
          if (firstLevel == refVar.getLevels()) {
            addLastLevel = true;
          }
        } else {
          // No extra level for 3D variables
          ASSERT(!extraLevel);
        }

        // Insert first level
        firstLevel_.insert({var.name(), firstLevel});
      }

      // Compute final number of levels and update first level index
      int refVarLevels = refVar.getLevels();
      if (addFirstLevel) {
        ++refVarLevels;
        for (auto & pair : firstLevel_) {
          ++pair.second;
        }
      }
      if (addLastLevel) {
        ++refVarLevels;
      }

      // Append group parameters
      groupInputVars_.push_back(groupVars);
      groupNames_.push_back(groupParams.groupName.value());
      groupInnerVars_.push_back(oops::Variables(
        {oops::Variable(groupNames_.back(), refVar.metaData(), refVarLevels)}));
      doCalibration_.push_back(groupParams.centralBlockParams().doCalibration());
      doRead_.push_back(groupParams.centralBlockParams().doRead());
      forceWrite_.push_back(groupParams.centralBlockParams().forceWrite.value());

      if (groupParams.outerBlocks.value()) {
        // Append group outer block chain
        oops::FieldSet4D fset4dXb(xb);
        oops::FieldSet4D fset4dFg(fg);
        groupOuterBlockChains_.push_back(std::make_unique<SaberOuterBlockChain>(outerGeom,
            groupInnerVars_.back(),
            fset4dXb,
            fset4dFg,
            covarConf,
            *groupParams.outerBlocks.value()));
      } else {
        // Append empty outer block chain
        std::unique_ptr<SaberOuterBlockChain> emptyOuterBlockChainPtr;
        groupOuterBlockChains_.push_back(std::move(emptyOuterBlockChainPtr));
      }

      // Set outer geometry data for central block
      const oops::GeometryData & currentOuterGeom = groupOuterBlockChains_.back() ?
                             groupOuterBlockChains_.back()->innerGeometryData() : outerGeom;

      // Set outer variables for central block
      const oops::Variables currentOuterVars = groupOuterBlockChains_.back() ?
                             groupOuterBlockChains_.back()->innerVars() : groupInnerVars_.back();

      // Append group central block
      groupCentralBlocks_.push_back(SaberCentralBlockFactory::create(currentOuterGeom,
                                                         currentOuterVars,
                                                         covarConf,
                                                         groupParams.centralBlockParams(),
                                                         xb,
                                                         fg));

      // Strategy-specific setup
      if (strategy_ == "duplicated and weighted") {
        // Prepare weights for the "duplicated and weighted" strategy

        // Allocation
        const size_t nv = groupVars.size();
        Eigen::MatrixXd wgt = Eigen::MatrixXd::Zero(nv, nv);

        // Set default weights
        const double defaultWeight = groupParams.defOffDiagWeight.value();
        for (size_t jvarJ = 0; jvarJ < groupVars.size(); ++jvarJ) {
          for (size_t jvarI = 0; jvarI < groupVars.size(); ++jvarI) {
            if (jvarJ == jvarI) {
              // Unit diagonal
              wgt(jvarJ, jvarI) = 1.0;
            } else {
              // Default off-diagonal weight
              wgt(jvarJ, jvarI) = defaultWeight;
            }
          }
        }

        // Set specific weights
        if (groupParams.offDiagWeights.value()) {
          for (const auto & specWeight : *groupParams.offDiagWeights.value()) {
            // Get variables pair and weight
            const std::vector<std::string> varPair = specWeight.varPair.value();
            ASSERT(varPair.size() == 2);
            const double weight = specWeight.weight.value();

            // Get variables pair indices
            const size_t jvarJ = groupVars.find(varPair[0]);
            const size_t jvarI = groupVars.find(varPair[1]);

            // Check that variables are different
            ASSERT(jvarJ != jvarI);

            // Set weight symmetrically
            wgt(jvarI, jvarJ) = weight;
            wgt(jvarJ, jvarI) = weight;
          }
        }

        // Cholesky decomposition
        wgtSqrt_.push_back(wgt.llt().matrixL());
        oops::Log::info() << "Info     : Weights square-root matrix for group "
                        << groupNames_.back() << " :" << std::endl;
        for (const auto row : wgtSqrt_.back().rowwise()) {
          oops::Log::info() << "Info     : " << row << std::endl;
        }
      }
    }
  }

  oops::Log::trace() << "SaberCentralBlock::SaberCentralBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::multiply(oops::FieldSet3D & fset3d) const {
  oops::Log::trace() << "SaberCentralBlock::multiply starting" << std::endl;

  if (strategy_ == "single") {
    // Single mode
    groupCentralBlocks_[0]->multiply(fset3d);
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate strategy or duplicated and weighted strategy
    if (strategy_ == "duplicated and weighted") {
      // Apply weights, adjoint
      applyWeightsAD(fset3d);
    }
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      for (const auto & var : groupInputVars_[igroup]) {
        // Create an empty FieldSet3D
        oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});

        // Clone input field
        auto inputField = fset3d[var.name()].clone();

        // Add field to empty FieldSet3D
        fset3dTmp.add(inputField);

        // Rename field with the name of the group
        inputField.rename(groupNames_[igroup]);

        // Multiply
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocksAD(fset3dTmp);
        groupCentralBlock->multiply(fset3dTmp);
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

        // Get output field
        const auto outputField = fset3dTmp[groupNames_[igroup]];

        // Copy content to input FieldSet3D
        const auto outputView = make_view<double, 2>(outputField);
        auto view = make_view<double, 2>(fset3d[var.name()]);
        view.assign(outputView);
      }
    }
    if (strategy_ == "duplicated and weighted") {
      // Apply weights
      applyWeights(fset3d);
    }
  } else if (strategy_ == "duplicated") {
    // Duplicated strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      // Create a FieldSet3D
      oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
      fset3dTmp.init(geometryData_.functionSpace(), groupInnerVars_[igroup]);
      fset3dTmp.zero();

      // Get input reference field
      auto inputRefField = fset3dTmp[groupNames_[igroup]];

      // Get input reference field view
      auto inputRefView = make_view<double, 2>(inputRefField);

      // Sum of fields
      for (const auto & var : groupInputVars_[igroup]) {
        // Get field
        const auto field = fset3d[var.name()];

        // Get field view
        const auto view = make_view<double, 2>(field);

        // Get first level
        const int firstLevel = firstLevel_.at(var.name());

        // Check whether the field is a 2D field added to a reference 3D field
        if ((inputRefField.levels() > 1) && (field.levels() == 1)) {
          // 2D level only
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            inputRefView(jnode, firstLevel) += view(jnode, 0);
          }
        } else {
          // All levels
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              inputRefView(jnode, firstLevel+jlevel) += view(jnode, jlevel);
            }
          }
        }
      }

      // Multiply
      if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocksAD(fset3dTmp);
      groupCentralBlock->multiply(fset3dTmp);
      if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

      // Get output reference field
      const auto outputRefField = fset3dTmp[groupNames_[igroup]];

      // Get output reference field view
      const auto outputRefView = make_view<double, 2>(outputRefField);

      // Split field
      for (const auto & var : groupInputVars_[igroup]) {
        // Get field
        auto field = fset3d[var.name()];

        // Get field view
        auto view = make_view<double, 2>(field);

        // Get first level
        const int firstLevel = firstLevel_.at(var.name());

        // Check whether the field is a 2D field added to a reference 3D field
        if ((outputRefField.levels() > 1) && (field.levels() == 1)) {
          // 2D level only
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            view(jnode, 0) = outputRefView(jnode, firstLevel);
          }
        } else {
          // All levels
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              view(jnode, jlevel) = outputRefView(jnode, firstLevel+jlevel);
            }
          }
        }
      }
    }
  } else if (strategy_ == "crossed") {
    // Crossed strategy

    // Initialization
    atlas::Field ctlVec = atlas::Field("genericCtlVec", make_datatype<double>(),
      make_shape(ctlVecSize()));

    // Apply multiplySqrtAD
    multiplySqrtAD(fset3d, ctlVec, 0);

    // Apply multiplySqrt
    multiplySqrt(ctlVec, fset3d, 0);
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  oops::Log::trace() << "SaberCentralBlock::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D SaberCentralBlock::variance() const {
  oops::Log::trace() << "SaberCentralBlock::variance starting" << std::endl;

  if (strategy_ == "single") {
    // Single mode: a single underlying central block carries the variance.
    return groupCentralBlocks_[0]->variance();
  } else {
    // Multivariate strategies are not yet supported because the variance
    // composition across groups depends on the strategy. Add coverage as
    // needed.
    throw eckit::NotImplemented("SaberCentralBlock::variance not implemented "
                                "for multivariate strategy '" + strategy_ + "'", Here());
  }
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::randomize(oops::FieldSet3D & fset3d) const {
  oops::Log::trace() << "SaberCentralBlock::randomize starting" << std::endl;

  if (strategy_ == "single") {
    // Single mode
    groupCentralBlocks_[0]->randomize(fset3d);
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate strategy or duplicated and weighted strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      for (const auto & var : groupInputVars_[igroup]) {
        // Create a FieldSet3D
        oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
        fset3dTmp.init(groupCentralBlock->geometryData().functionSpace(),
                       groupCentralBlock->centralVars());

        // Randomize
        groupCentralBlock->randomize(fset3dTmp);
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

        // Get output field
        const auto outputField = fset3dTmp[groupNames_[igroup]];

        // Copy content to input FieldSet3D
        const auto outputView = make_view<double, 2>(outputField);
        auto view = make_view<double, 2>(fset3d[var.name()]);
        view.assign(outputView);
      }
    }
    if (strategy_ == "duplicated and weighted") {
      // Apply weights
      applyWeights(fset3d);
    }
  } else if ((strategy_ == "duplicated") || (strategy_ == "crossed")) {
    // Duplicated strategy or crossed strategy
    atlas::Field cv;
    if (strategy_ == "crossed") {
      // Create control vector
      cv = atlas::Field("genericCtlVec", atlas::array::make_datatype<double>(),
        atlas::array::make_shape(ctlVecSize()));

      // Generate random control vector
      randomCtlVec(cv, 0);
    }
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      // Create a FieldSet3D
      oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
      fset3dTmp.init(groupCentralBlock->geometryData().functionSpace(),
                     groupCentralBlock->centralVars());

      // Randomize
      if (strategy_ == "crossed") {
        // Square-root multiply
        groupCentralBlock->multiplySqrt(cv, fset3d, 0);
      } else if (strategy_ == "duplicated") {
        // Independent randomization for each group
        groupCentralBlock->randomize(fset3dTmp);
      }
      if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

      // Get output reference field
      const auto outputRefField = fset3dTmp[groupNames_[igroup]];

      // Get output reference field view
      const auto outputRefView = make_view<double, 2>(outputRefField);

      // Split field
      for (const auto & var : groupInputVars_[igroup]) {
        // Get field
        auto field = fset3d[var.name()];

        // Get field view
        auto view = make_view<double, 2>(field);

        // Get first level
        const int firstLevel = firstLevel_.at(var.name());

        // Check whether the field is a 2D field added to a reference 3D field
        if ((outputRefField.levels() > 1) && (field.levels() == 1)) {
          // 2D level only
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            view(jnode, 0) = outputRefView(jnode, firstLevel);
          }
        } else {
          // All levels
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              view(jnode, jlevel) = outputRefView(jnode, firstLevel+jlevel);
            }
          }
        }
      }
    }
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  oops::Log::trace() << "SaberCentralBlock::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::read() {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    if (doRead_[igroup]) {
      groupCentralBlocks_[igroup]->read();
    }
  }
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::directCalibration(const oops::FieldSets & fsets) {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    if (doCalibration_[igroup]) {
      if (groupOuterBlockChains_[igroup]) {
        // Applying the left inverse of the outer block chain would require a deep-copy
        // of the ensemble, not allowed here. If needed, ensemble-based calibration should be done
        // outside of the multivariate central block framework.
        ASSERT(fsets.size() == 0);
      }
      groupCentralBlocks_[igroup]->directCalibration(fsets);
    }
  }
}


// -----------------------------------------------------------------------------

void SaberCentralBlock::iterativeCalibrationInit() {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    if (doCalibration_[igroup]) {
      ASSERT(!groupOuterBlockChains_[igroup]);
      groupCentralBlocks_[igroup]->iterativeCalibrationInit();
    }
  }
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    if (doCalibration_[igroup]) {
      ASSERT(!groupOuterBlockChains_[igroup]);
      groupCentralBlocks_[igroup]->iterativeCalibrationUpdate(fset);
    }
  }
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::iterativeCalibrationFinal() {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    if (doCalibration_[igroup]) {
      ASSERT(!groupOuterBlockChains_[igroup]);
      groupCentralBlocks_[igroup]->iterativeCalibrationFinal();
    }
  }
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::write() const {
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    groupCentralBlocks_[igroup]->write();
  }
}

// -----------------------------------------------------------------------------

size_t SaberCentralBlock::ctlVecSize() const {
  // Initialize control vector size
  size_t ctlVecSize = 0;

  if (strategy_ == "single") {
    // Single mode
    ctlVecSize += groupCentralBlocks_[0]->ctlVecSize();
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Add the group control vector size for each variable
      ctlVecSize += groupCentralBlocks_[igroup]->ctlVecSize()*groupInputVars_[igroup].size();
    }
  } else if (strategy_ == "duplicated") {
    // Duplicated strategy
    for (const auto & groupCentralBlock : groupCentralBlocks_) {
      // Add the group control vector
      ctlVecSize += groupCentralBlock->ctlVecSize();
    }
  } else if (strategy_ == "crossed") {
    // Crossed strategy
    ctlVecSize += groupCentralBlocks_[0]->ctlVecSize();
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  // Return control vector size
  return ctlVecSize;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::randomCtlVec(atlas::Field & cv,
                                     const size_t & offset) const {
  oops::Log::trace() << "SaberCentralBlock::randomCtlVec starting" << std::endl;

  // Initialize index
  size_t index = offset;

  if (strategy_ == "single") {
    // Single mode
    groupCentralBlocks_[0]->randomCtlVec(cv, offset);
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate or duplicated and weighted strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      for (size_t jvar = 0; jvar < groupInputVars_[igroup].size(); ++jvar) {
        // Get random control vector
        groupCentralBlocks_[igroup]->randomCtlVec(cv, index);

        // Update index
        index += groupCentralBlocks_[igroup]->ctlVecSize();
      }
    }
  } else if (strategy_ == "duplicated") {
    // Duplicated strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get random control vector
      groupCentralBlocks_[igroup]->randomCtlVec(cv, index);

      // Update index
      index += groupCentralBlocks_[igroup]->ctlVecSize();
    }
  } else if (strategy_ == "crossed") {
    // Crossed strategy
    groupCentralBlocks_[0]->randomCtlVec(cv, 0);
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  oops::Log::trace() << "SaberCentralBlock::randomCtlVec done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::multiplySqrt(const atlas::Field & cv,
                                     oops::FieldSet3D & fset3d,
                                     const size_t & offset) const {
  oops::Log::trace() << "SaberCentralBlock::multiplySqrt starting" << std::endl;

  // Initialize index
  size_t index = offset;

  if (strategy_ == "single") {
    // Single mode
    groupCentralBlocks_[0]->multiplySqrt(cv, fset3d, offset);
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate strategy or duplicated and weighted strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      for (const auto & var : groupInputVars_[igroup]) {
        // Create a FieldSet3D
        oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
        fset3dTmp.init(groupCentralBlock->geometryData().functionSpace(),
                       groupCentralBlock->centralVars());

        // Square-root multiply
        groupCentralBlock->multiplySqrt(cv, fset3dTmp, index);
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

        // Update index
        index += groupCentralBlock->ctlVecSize();

        // Get output field
        const auto outputField = fset3dTmp[groupNames_[igroup]];

        // Copy content to input FieldSet3D
        const auto outputView = make_view<double, 2>(outputField);
        auto view = make_view<double, 2>(fset3d[var.name()]);
        view.assign(outputView);
      }
    }
    if (strategy_ == "duplicated and weighted") {
      // Apply weights
      applyWeights(fset3d);
    }
  } else if ((strategy_ == "duplicated") || (strategy_ == "crossed")) {
    // Duplicated strategy or crossed strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      // Create a FieldSet3D
      oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
      fset3dTmp.init(groupCentralBlock->geometryData().functionSpace(),
                     groupCentralBlock->centralVars());

      // Square-root multiply
      groupCentralBlock->multiplySqrt(cv, fset3dTmp, index);
      if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocks(fset3dTmp);

      if (strategy_ == "duplicated") {
        // Update index
        index += groupCentralBlock->ctlVecSize();
      }

      // Get output reference field
      const auto outputRefField = fset3dTmp[groupNames_[igroup]];

      // Get output reference field view
      const auto outputRefView = make_view<double, 2>(outputRefField);

      // Split field
      for (const auto & var : groupInputVars_[igroup]) {
        // Get field
        auto field = fset3d[var.name()];

        // Get field view
        auto view = make_view<double, 2>(field);

        // Get first level
        const int firstLevel = firstLevel_.at(var.name());

        // Check whether the field is a 2D field added to a reference 3D field
        if ((outputRefField.levels() > 1) && (field.levels() == 1)) {
          // 2D level only
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            view(jnode, 0) = outputRefView(jnode, firstLevel);
          }
        } else {
          // All levels
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              view(jnode, jlevel) = outputRefView(jnode, firstLevel+jlevel);
            }
          }
        }
      }
    }
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  oops::Log::trace() << "SaberCentralBlock::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::multiplySqrtAD(const oops::FieldSet3D & fset3d,
                                       atlas::Field & cv,
                                       const size_t & offset) const {
  oops::Log::trace() << "SaberCentralBlock::multiplySqrtAD starting" << std::endl;

  // Initialize index
  size_t index = offset;

  // Initialize control variable
  auto ctlVecView = make_view<double, 1>(cv);
  for (size_t jnode = 0; jnode < ctlVecSize(); ++jnode) {
    ctlVecView(index+jnode) = 0.0;
  }

  if (strategy_ == "single") {
    // Single mode
    groupCentralBlocks_[0]->multiplySqrtAD(fset3d, cv, offset);
  } else if ((strategy_ == "univariate") || (strategy_ == "duplicated and weighted")) {
    // Univariate strategy or duplicated and weighted strategy
    oops::FieldSet3D fset3dCopy(fset3d);
    if (strategy_ == "duplicated and weighted") {
      // Apply weights, adjoint
      applyWeightsAD(fset3dCopy);
    }
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      for (const auto & var : groupInputVars_[igroup]) {
        // Create an empty FieldSet3D
        oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});

        // Clone field
        auto field = fset3dCopy[var.name()].clone();

        // Rename field with the name of the group
        field.rename(groupNames_[igroup]);

        // Add field
        fset3dTmp.add(field);

        // Square-root adjoint multiply
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocksAD(fset3dTmp);
        groupCentralBlock->multiplySqrtAD(fset3dTmp, cv, index);

        // Update index
        index += groupCentralBlock->ctlVecSize();
      }
    }
  } else if ((strategy_ == "duplicated") || (strategy_ == "crossed")) {
    // Duplicated strategy or crossed strategy
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      // Get outer block chain
      const auto & groupOuterBlockChain = groupOuterBlockChains_[igroup];

      // Get central block
      const auto & groupCentralBlock = groupCentralBlocks_[igroup];

      // Create a FieldSet3D
      oops::FieldSet3D fset3dTmp({fset3d.validTime(), fset3d.commGeom()});
      fset3dTmp.init(geometryData_.functionSpace(), groupInnerVars_[igroup]);

      // Get input reference field
      auto inputRefField = fset3dTmp[groupNames_[igroup]];

      // Get reference field view
      auto inputRefView = make_view<double, 2>(inputRefField);
      inputRefView.assign(0.0);

      // Sum of fields
      for (const auto & var : groupInputVars_[igroup]) {
        // Get field
        const auto field = fset3d[var.name()];

        // Get field view
        const auto view = make_view<double, 2>(field);

        // Get first level
        const int firstLevel = firstLevel_.at(var.name());

        // Check whether the field is a 2D field added to a reference 3D field
        if ((inputRefField.levels() > 1) && (field.levels() == 1)) {
          // 2D level only
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            inputRefView(jnode, firstLevel) += view(jnode, 0);
          }
        } else {
          // All levels
          for (int jnode = 0; jnode < field.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
              inputRefView(jnode, firstLevel+jlevel) += view(jnode, jlevel);
            }
          }
        }
      }

      if (strategy_ == "duplicated") {
        // Square-root adjoint multiply
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocksAD(fset3dTmp);
        groupCentralBlock->multiplySqrtAD(fset3dTmp, cv, index);

        // Update index
        index += groupCentralBlock->ctlVecSize();
      } else if (strategy_ == "crossed") {
        // Create temporary control vector
        atlas::Field ctlVecTmp = atlas::Field("genericCtlVec", make_datatype<double>(),
          make_shape(groupCentralBlock->ctlVecSize()));

        // Square-root adjoint multiply
        if (groupOuterBlockChain) groupOuterBlockChain->applyOuterBlocksAD(fset3dTmp);
        groupCentralBlock->multiplySqrtAD(fset3dTmp, ctlVecTmp, 0);

        // Add control vector contribution
        const auto ctlVecTmpView = make_view<double, 1>(ctlVecTmp);
        for (int jnode = 0; jnode < ctlVecTmp.shape(0); ++jnode) {
          ctlVecView(index+jnode) += ctlVecTmpView(jnode);
        }
      }
    }
  } else {
    throw eckit::Exception("invalid multivariate strategy", Here());
  }

  oops::Log::trace() << "SaberCentralBlock::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::calibrateBlock(const oops::FieldSet4D & fset4dXb) {
  oops::Log::trace() << "SaberCentralBlock::calibrateBlock starting" << std::endl;

  // Create empty ensemble
  std::vector<util::DateTime> dates;
  std::vector<int> ensmems;
  oops::FieldSets fsetEns(dates, fset4dXb.commTime(), ensmems, fset4dXb.commEns());

  // Direct calibration
  oops::Log::info() << "Info     : Direct calibration (without ensemble)" << std::endl;
  this->directCalibration(fsetEns);

  oops::Log::trace() << "SaberCentralBlock::calibrateBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::adjointTest(const double & globalAdjointTolerance) const {
  oops::Log::trace() << "SaberCentralBlock::adjointTest starting" << std::endl;

  // Separate test for each group
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    // Override adjoint tolerance if specified in the configuration
    double adjointTolerance = globalAdjointTolerance;
    if (params_.singleBlock.value() && params_.singleBlock.value()->adjointTolerance.value()) {
      adjointTolerance = *params_.singleBlock.value()->adjointTolerance.value();
    } else if (params_.groups.value()) {
      const auto & groupParams = params_.groups.value().get()[igroup];
      if (groupParams.centralBlockParams().adjointTolerance.value()) {
        adjointTolerance = *groupParams.centralBlockParams().adjointTolerance.value();
      }
    }

    // Get geometry data
    const auto & geomData = groupCentralBlocks_[igroup]->geometryData();

    // Get central variables
    const auto & centralVars = groupCentralBlocks_[igroup]->centralVars();

    // Create random FieldSets
    oops::FieldSet3D fset1 = oops::randomFieldSet3D(validTime_,
                                                    geomData.comm(),
                                                    geomData.functionSpace(),
                                                    centralVars);
    oops::FieldSet3D fset2 = oops::randomFieldSet3D(validTime_,
                                                    geomData.comm(),
                                                    geomData.functionSpace(),
                                                    centralVars);

    // Copy FieldSets
    oops::FieldSet3D fset1Save(fset1);
    oops::FieldSet3D fset2Save(fset2);

    // Apply forward multiplication only (self-adjointness test)
    groupCentralBlocks_[igroup]->multiply(fset1);
    groupCentralBlocks_[igroup]->multiply(fset2);

    // Compute adjoint test
    const double dp1 = fset1.dot_product_with(fset2Save, centralVars);
    const double dp2 = fset2.dot_product_with(fset1Save, centralVars);
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: (Ax)^t y = " << dp1
                      << ": x^t (Ay) = " << dp2 << " : adjoint tolerance = "
                      << adjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for block " << groupCentralBlocks_[igroup]->blockName();
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < adjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Adjoint test failure for block " +
        groupCentralBlocks_[igroup]->blockName(), Here());
    }
  }

  // Test for all the groups at once
  if (params_.groups.value()) {
    // Override adjoint tolerance if specified in the configuration
    double adjointTolerance = globalAdjointTolerance;
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      const auto & groupParams = params_.groups.value().get()[igroup];
      if (groupParams.centralBlockParams().adjointTolerance.value()) {
        adjointTolerance = std::max(adjointTolerance,
          *groupParams.centralBlockParams().adjointTolerance.value());
      }
    }

    // Create random FieldSets
    oops::FieldSet3D fset1 = oops::randomFieldSet3D(validTime_,
                                                    geometryData_.comm(),
                                                    geometryData_.functionSpace(),
                                                    outerVars_);
    oops::FieldSet3D fset2 = oops::randomFieldSet3D(validTime_,
                                                    geometryData_.comm(),
                                                    geometryData_.functionSpace(),
                                                    outerVars_);

    // Copy FieldSets
    oops::FieldSet3D fset1Save(fset1);
    oops::FieldSet3D fset2Save(fset2);

    // Apply forward multiplication only (self-adjointness test)
    this->multiply(fset1);
    this->multiply(fset2);

    // Compute adjoint test
    const double dp1 = fset1.dot_product_with(fset2Save, outerVars_);
    const double dp2 = fset2.dot_product_with(fset1Save, outerVars_);
    oops::Log::info() << std::setprecision(16) << "Info     : Adjoint test: (Ax)^t y = " << dp1
                      << ": x^t (Ay) = " << dp2 << " : adjoint tolerance = "
                      << adjointTolerance << std::endl;
    oops::Log::test() << "Adjoint test for the multivariate central block";
    if (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < adjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Adjoint test failure for the multivariate central block", Here());
    }
  }

  oops::Log::trace() << "SaberCentralBlock::adjointTest done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::sqrtTest(const double & globalSqrtTolerance) const {
  oops::Log::trace() << "SaberOuterBlockBase::sqrtTest starting" << std::endl;

  // Separate test for each group
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    // Override square root tolerance if specified in the configuration
    double sqrtTolerance = globalSqrtTolerance;
    if (params_.singleBlock.value() && params_.singleBlock.value()->sqrtTolerance.value()) {
      sqrtTolerance = *params_.singleBlock.value()->sqrtTolerance.value();
    } else if (params_.groups.value()) {
      const auto & groupParams = params_.groups.value().get()[igroup];
      if (groupParams.centralBlockParams().sqrtTolerance.value()) {
        sqrtTolerance = *groupParams.centralBlockParams().sqrtTolerance.value();
      }
    }

    // Get geometry data
    const auto & geomData = groupCentralBlocks_[igroup]->geometryData();

    // Get central variables
    const auto & centralVars = groupCentralBlocks_[igroup]->centralVars();

    // Create FieldSet
    oops::FieldSet3D fset = oops::randomFieldSet3D(validTime_,
                                                   geomData.comm(),
                                                   geomData.functionSpace(),
                                                   centralVars);

    // Copy FieldSet
    oops::FieldSet3D fsetSave(fset);

    // Create control vector for this group
    const size_t ctlVecSize = groupCentralBlocks_[igroup]->ctlVecSize();
    oops::Log::info() << "Info     : Control vector size for block "
                      << groupCentralBlocks_[igroup]->blockName() << ": "
                      << ctlVecSize << std::endl;
    atlas::Field ctlVec = atlas::Field("genericCtlVec",
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::make_shape(ctlVecSize));
    size_t seed = 7;  // To avoid impact on future random generator calls
    util::NormalDistribution<double> dist(ctlVecSize, 0.0, 1.0, seed);
    auto view = atlas::array::make_view<double, 1>(ctlVec);
    for (size_t jnode = 0; jnode < ctlVecSize; ++jnode) {
      view(jnode) = dist[jnode];
    }

    // Copy control vector
    atlas::Field ctlVecSave = atlas::Field("genericCtlVec",
                                           atlas::array::make_datatype<double>(),
                                           atlas::array::make_shape(ctlVecSize));
    auto viewSave = atlas::array::make_view<double, 1>(ctlVecSave);
    viewSave.assign(view);

    // Apply square-root multiplication
    groupCentralBlocks_[igroup]->multiplySqrt(ctlVecSave, fset, 0);

    // Apply square-root adjoint multiplication
    groupCentralBlocks_[igroup]->multiplySqrtAD(fsetSave, ctlVec, 0);

    // Compute adjoint test
    const double dp1 = fset.dot_product_with(fsetSave, centralVars);
    double dp2 = 0.0;
    for (size_t jnode = 0; jnode < ctlVecSize; ++jnode) {
      dp2 += view(jnode)*viewSave(jnode);
    }
    geomData.comm().allReduceInPlace(dp2, eckit::mpi::sum());
    oops::Log::info() << std::setprecision(16) << "Info     : Square-root test: y^t (Ux) = " << dp1
                      << ": x^t (U^t y) = " << dp2 << " : square-root tolerance = "
                      << sqrtTolerance << std::endl;
    const bool adjComparison = (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < sqrtTolerance);

    // Apply square-root multiplication
    groupCentralBlocks_[igroup]->multiplySqrt(ctlVec, fset, 0);

    // Apply full multiplication
    groupCentralBlocks_[igroup]->multiply(fsetSave);

    // Check that the fieldsets are similar within tolerance
    const bool sqrtComparison = fset.compare_with(fsetSave, sqrtTolerance,
                                                  util::ToleranceType::relative);
    if (sqrtComparison) {
      oops::Log::info() << "Info     : Square-root test passed: U U^t x == B x" << std::endl;
    } else {
      oops::Log::info() << "Info     : Square-root test failed: U U^t x != B x" << std::endl;
    }

    // Print results
    oops::Log::test() << "Square-root test for block " << groupCentralBlocks_[igroup]->blockName();
    if (adjComparison && sqrtComparison) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Square-root test failure for block "
         + groupCentralBlocks_[igroup]->blockName(), Here());
    }
  }

  // Test for all the groups at once
  if (params_.groups.value()) {
    // Override square root tolerance if specified in the configuration
    double sqrtTolerance = globalSqrtTolerance;
    for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
      const auto & groupParams = params_.groups.value().get()[igroup];
      if (groupParams.centralBlockParams().sqrtTolerance.value()) {
        sqrtTolerance = std::max(sqrtTolerance,
          *groupParams.centralBlockParams().sqrtTolerance.value());
      }
    }

    // Create FieldSet
    oops::FieldSet3D fset = oops::randomFieldSet3D(validTime_,
                                                   geometryData_.comm(),
                                                   geometryData_.functionSpace(),
                                                   outerVars_);

    // Copy FieldSet
    oops::FieldSet3D fsetSave(fset);

    // Create control vector for this group
    const size_t ctlVecSize = this->ctlVecSize();
    oops::Log::info() << "Info     : Control vector size for the multivariate central block: "
                      << ctlVecSize << std::endl;
    atlas::Field ctlVec = atlas::Field("genericCtlVec",
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::make_shape(ctlVecSize));
    size_t seed = 7;  // To avoid impact on future random generator calls
    util::NormalDistribution<double> dist(ctlVecSize, 0.0, 1.0, seed);
    auto view = atlas::array::make_view<double, 1>(ctlVec);
    for (size_t jnode = 0; jnode < ctlVecSize; ++jnode) {
      view(jnode) = dist[jnode];
    }

    // Copy control vector
    atlas::Field ctlVecSave = atlas::Field("genericCtlVec",
                                           atlas::array::make_datatype<double>(),
                                           atlas::array::make_shape(ctlVecSize));
    auto viewSave = atlas::array::make_view<double, 1>(ctlVecSave);
    viewSave.assign(view);

    // Apply square-root multiplication
    this->multiplySqrt(ctlVecSave, fset, 0);

    // Apply square-root adjoint multiplication
    this->multiplySqrtAD(fsetSave, ctlVec, 0);

    // Compute adjoint test
    const double dp1 = fset.dot_product_with(fsetSave, outerVars_);
    double dp2 = 0.0;
    for (size_t jnode = 0; jnode < ctlVecSize; ++jnode) {
      dp2 += view(jnode)*viewSave(jnode);
    }
    geometryData_.comm().allReduceInPlace(dp2, eckit::mpi::sum());
    oops::Log::info() << std::setprecision(16) << "Info     : Square-root test: y^t (Ux) = " << dp1
                      << ": x^t (U^t y) = " << dp2 << " : square-root tolerance = "
                      << sqrtTolerance << std::endl;
    const bool adjComparison = (std::abs(dp1-dp2)/std::abs(0.5*(dp1+dp2)) < sqrtTolerance);

    // Apply square-root multiplication
    this->multiplySqrt(ctlVec, fset, 0);

    // Apply full multiplication
    this->multiply(fsetSave);

    // Check that the fieldsets are similar within tolerance
    const bool sqrtComparison = fset.compare_with(fsetSave, sqrtTolerance,
                                                  util::ToleranceType::relative);
    if (sqrtComparison) {
      oops::Log::info() << "Info     : Square-root test passed: U U^t x == B x" << std::endl;
    } else {
      oops::Log::info() << "Info     : Square-root test failed: U U^t x != B x" << std::endl;
    }

    // Print results
    oops::Log::test() << "Square-root test for the multivariate central block";
    if (adjComparison && sqrtComparison) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      throw eckit::Exception("Square-root test failure for the multivariate central block", Here());
    }
  }
  oops::Log::trace() << "SaberCentralBlock::sqrtTest done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::applyWeights(oops::FieldSet3D & fset3d) const {
  oops::Log::trace() << "SaberCentralBlock::applyWeights starting" << std::endl;

  // Check strategy
  ASSERT(strategy_ == "duplicated and weighted");

  // Deep copy
  oops::FieldSet3D fset3dCopy(fset3d);

  // Set input to zero
  fset3d.zero();

  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    for (size_t jvarI = 0; jvarI < groupInputVars_[igroup].size(); ++jvarI) {
      // Get variable I
      const auto varI = groupInputVars_[igroup][jvarI];

      // Get field I
      auto fieldI = fset3d[varI.name()];

      // Get field I view
      auto viewI = make_view<double, 2>(fieldI);

      for (size_t jvarJ = 0; jvarJ <= jvarI; ++jvarJ) {
        // Get variable J
        const auto varJ = groupInputVars_[igroup][jvarJ];

        // Get field J
        const auto fieldJ = fset3dCopy[varJ.name()];

        // Get field J view
        const auto viewJ = make_view<double, 2>(fieldJ);

        // Sum weighted off-diagonal fields
        for (int jnode = 0; jnode < fieldJ.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < fieldJ.shape(1); ++jlevel) {
            viewI(jnode, jlevel) += wgtSqrt_[igroup](jvarI, jvarJ)*viewJ(jnode, jlevel);
          }
        }
      }
    }
  }

  oops::Log::trace() << "SaberCentralBlock::applyWeights done" << std::endl;
}

// -----------------------------------------------------------------------------

void SaberCentralBlock::applyWeightsAD(oops::FieldSet3D & fset3d) const {
  oops::Log::trace() << "SaberCentralBlock::applyWeightsAD starting" << std::endl;

  // Check strategy
  ASSERT(strategy_ == "duplicated and weighted");

  // Deep copy
  oops::FieldSet3D fset3dCopy(fset3d);

  // Set input to zero
  fset3d.zero();

  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    for (size_t jvarI = 0; jvarI < groupInputVars_[igroup].size(); ++jvarI) {
      // Get variable I
      const auto varI = groupInputVars_[igroup][jvarI];

      // Get field I
      auto fieldI = fset3d[varI.name()];

      // Get field I view
      auto viewI = make_view<double, 2>(fieldI);

      for (size_t jvarJ = jvarI; jvarJ < groupInputVars_[igroup].size(); ++jvarJ) {
        // Get variable J
        const auto varJ = groupInputVars_[igroup][jvarJ];

        // Get field J
        const auto fieldJ = fset3dCopy[varJ.name()];

        // Get field J view
        const auto viewJ = make_view<double, 2>(fieldJ);

        // Sum weighted off-diagonal fields
        for (int jnode = 0; jnode < fieldJ.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < fieldJ.shape(1); ++jlevel) {
            viewI(jnode, jlevel) += wgtSqrt_[igroup](jvarJ, jvarI)*viewJ(jnode, jlevel);
          }
        }
      }
    }
  }

  oops::Log::trace() << "SaberCentralBlock::applyWeightsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

