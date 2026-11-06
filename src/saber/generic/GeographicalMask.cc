/*
 * (C) Copyright 2026 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <omp.h>

#include <algorithm>
#include <string>

#include "saber/generic/GeographicalMask.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

#include "saber/blocks/OuterBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static OuterBlockMaker<GeographicalMask> makerGeographicalMask_("GeographicalMask");

// -----------------------------------------------------------------------------

GeographicalMask::GeographicalMask(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : OuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerGeometryData_(outerGeometryData),
    activeOuterVars_(params.activeVars.value().get_value_or(outerVars)),
    params_(params.calibration.value() ? *params.calibration.value() : *params.read.value()) {
  oops::Log::trace() << classname() << "::GeographicalMask starting" << std::endl;

  // Define inner variables
  for (const auto & outerVar : outerVars) {
    if (activeOuterVars_.has(outerVar)) {
      // Add active inner variable
      for (const auto & maskParams : params_.masksParams.value()) {
        const oops::Variable innerVar(outerVar.name() + maskParams.suffix.value(),
          outerVar.metaData(), outerVar.getLevels());
        ASSERT(!innerVars_.has(innerVar.name()));
        innerVars_.push_back(innerVar);
      }
    } else {
      // Add inactive inner variables
      innerVars_.push_back(outerVars[outerVar.name()]);
    }
  }

  // Check consistency
  if (params_.inputMaskAtlasFileConf.value() && params_.inputMaskModelFileConf.value()) {
    throw eckit::UserError("cannot read mask with both ATLAS and model readers", Here());
  }
  if (params_.inputWgtAtlasFileConf.value() && params_.inputWgtModelFileConf.value()) {
    throw eckit::UserError("cannot read weight with both ATLAS and model readers", Here());
  }
  if ((params_.inputMaskAtlasFileConf.value() || params_.inputMaskModelFileConf.value())
    && (params_.inputWgtAtlasFileConf.value() || params_.inputWgtModelFileConf.value())) {
    throw eckit::UserError("cannot read both mask and weight at the same time", Here());
  }

  if (params_.filterParams.value()) {
    // Save xb and fg for the filter constructor
    fset4dXb_ = std::make_unique<oops::FieldSet4D>(xb);
    fset4dFg_ = std::make_unique<oops::FieldSet4D>(fg);
  }

  if (params.calibration.value()) {
    // Create mask fieldset
    maskFset_ = std::make_unique<oops::FieldSet3D>(xb.validTime(), innerGeometryData_.comm());
  }

  // Create weight fieldset
  wgtFset_ = std::make_unique<oops::FieldSet3D>(xb.validTime(), innerGeometryData_.comm());

  oops::Log::trace() << classname() << "::GeographicalMask done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create fieldset
  atlas::FieldSet outerFset;

  for (const auto & outerVar : outerVars()) {
    if (activeOuterVars_.has(outerVar)) {
      // Create outer field
      atlas::Field outerField = outerGeometryData().functionSpace().createField<double>(
        atlas::option::name(outerVar.name()) | atlas::option::levels(outerVar.getLevels())
        | atlas::option::halo(1));

      // Get and initialize outer view
      auto outerView = atlas::array::make_view<double, 2>(outerField);
      outerView.assign(0.0);

      for (const auto & maskParams : params_.masksParams.value()) {
        // Get inner view
        const auto innerView = atlas::array::make_view<double, 2>(
          fset[outerVar.name() + maskParams.suffix.value()]);

        // Get weight field
        const auto wgtField = (*wgtFset_)[maskParams.suffix.value()];

        // Weight field is either 3D with the same number as outer variables, or 2D with one level
        ASSERT((wgtField.levels() == outerField.levels()) || (wgtField.levels() == 1));

        // Get weight view
        const auto wgtView = atlas::array::make_view<double, 2>(wgtField);

        // Add with a weight
        for (int jnode = 0; jnode < outerField.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < outerField.shape(1); ++jlevel) {
            outerView(jnode, jlevel) += innerView(jnode, jlevel)
              *wgtView(jnode, std::min(jlevel, wgtField.levels()-1));
          }
        }
      }

      // Add active inner variables
      outerFset.add(outerField);
    } else {
      // Add inactive inner variables
      outerFset.add(fset[outerVar.name()]);
    }
  }

  // Copy outer fieldset
  fset.fieldSet() = outerFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // Create fieldset
  atlas::FieldSet innerFset;

  for (const auto & outerVar : outerVars()) {
    if (activeOuterVars_.has(outerVar)) {
      // Get outer view
      const auto outerView = atlas::array::make_view<double, 2>(fset[outerVar.name()]);

      for (const auto & maskParams : params_.masksParams.value()) {
        // Create inner field
        atlas::Field innerField = innerGeometryData_.functionSpace().createField<double>(
          atlas::option::name(outerVar.name() + maskParams.suffix.value())
          | atlas::option::levels(outerVar.getLevels()) | atlas::option::halo(1));

        // Get inner view
        auto innerView = atlas::array::make_view<double, 2>(innerField);

        // Get weight field
        const auto wgtField = (*wgtFset_)[maskParams.suffix.value()];

        // Weight field is either 3D with the same number as outer variables, or 2D with one level
        ASSERT((wgtField.levels() == innerField.levels()) || (wgtField.levels() == 1));

        // Get weight view
        const auto wgtView = atlas::array::make_view<double, 2>(wgtField);

        // Copy and apply weight
        for (int jnode = 0; jnode < innerField.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < innerField.shape(1); ++jlevel) {
            innerView(jnode, jlevel) = outerView(jnode, jlevel)
              *wgtView(jnode, std::min(jlevel, wgtField.levels()-1));
          }
        }

        // Add inactive inner variables
        innerFset.add(innerField);
      }
    } else {
      // Add inactive inner variables
      innerFset.add(fset[outerVar.name()]);
    }
  }

  // Copy inner fieldset
  fset.fieldSet() = innerFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> GeographicalMask::getReadConfs()
  const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (params_.inputMaskModelFileConf.value()) {
    inputs.push_back(std::make_pair("mask", *params_.inputMaskModelFileConf.value()));
  }
  if (params_.inputWgtModelFileConf.value()) {
    inputs.push_back(std::make_pair("weight", *params_.inputWgtModelFileConf.value()));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void GeographicalMask::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  if (params_.inputMaskModelFileConf.value()) {
    // Get mask from input file
    ASSERT(fsetVec.size() == 1);
    maskFset_->deepCopy(fsetVec[0]);
  }
  if (params_.inputWgtModelFileConf.value()) {
    // Get weight from input file
    ASSERT(fsetVec.size() == 1);
    wgtFset_->deepCopy(fsetVec[0]);
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::directCalibration(const oops::FieldSets &) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  if (params_.inputMaskAtlasFileConf.value()) {
    // Read ATLAS mask file
    ASSERT(params_.maskVariable.value());
    ASSERT(params_.maskLevels.value());
    eckit::LocalConfiguration varConf;
    varConf.set("levels", *params_.maskLevels.value());
    const oops::Variables maskVars({oops::Variable(*params_.maskVariable.value(), varConf)});
    maskFset_->read(innerGeometryData_.functionSpace(),
                    maskVars,
                    *params_.inputMaskAtlasFileConf.value());
  }

  // Get ghost view
  const auto ghostView = atlas::array::make_view<int, 1>(
    outerGeometryData().functionSpace().ghost());

  // Get mask field
  ASSERT(params_.maskVariable.value());
  const auto & maskField = (*maskFset_)[*params_.maskVariable.value()];
  oops::Log::test() << "Mask norm: "
    << util::normField(maskField, outerGeometryData().comm()) << std::endl;

  // Get mask view
  const auto & maskView = atlas::array::make_view<double, 2>(maskField);

  for (const auto & maskParams : params_.masksParams.value()) {
    // Get relational operator and
    ASSERT(maskParams.relationalOperator.value());
    ASSERT(maskParams.threshold.value());
    const RelationalOperator & relationalOperator = *maskParams.relationalOperator.value();
    const double & threshold = *maskParams.threshold.value();

    // Create weight field
    atlas::Field wgtField = outerGeometryData().functionSpace().createField<double>(
      atlas::option::name(maskParams.suffix.value())
      | atlas::option::levels(maskField.levels()) | atlas::option::halo(1));
    auto wgtView = atlas::array::make_view<double, 2>(wgtField);

    // Check dimensions consistency
    ASSERT(maskField.shape() == wgtField.shape());

    // Fill weight field
    for (int jlevel = 0; jlevel < wgtField.shape(1); ++jlevel) {
      // Initialize counters
      int countInMask = 0;
      int countNotInMask = 0;

      for (int jnode = 0; jnode < wgtField.shape(0); ++jnode) {
        // Is point in mask?
        bool inMask = false;
        switch (relationalOperator) {
          case RelationalOperator::LESS_THAN: {
            inMask = (maskView(jnode, jlevel) < threshold);
            break;
          }
          case RelationalOperator::LESS_THAN_OR_EQUAL: {
            inMask = (maskView(jnode, jlevel) <= threshold);
            break;
          }
          case RelationalOperator::EQUAL: {
            inMask = (maskView(jnode, jlevel) == threshold);
            break;
          }
          case RelationalOperator::NOT_EQUAL: {
            inMask = (maskView(jnode, jlevel) != threshold);
            break;
          }
          case RelationalOperator::GREATER_THAN_OR_EQUAL: {
            inMask = (maskView(jnode, jlevel) >= threshold);
            break;
          }
          case RelationalOperator::GREATER_THAN: {
            inMask = (maskView(jnode, jlevel) > threshold);
            break;
          }
          default: {
            throw eckit::Exception("directCalibration: wrong relational operator", Here());
          }
        }

        // Set weight
        if (inMask) {
          wgtView(jnode, jlevel) = 1.0;
          if (ghostView(jnode) == 0) {
            ++countInMask;
          }
        } else {
          wgtView(jnode, jlevel) = 0.0;
          if (ghostView(jnode) == 0) {
            ++countNotInMask;
          }
        }
      }

      // Gather counter
      outerGeometryData().comm().allReduceInPlace(countInMask, eckit::mpi::sum());
      outerGeometryData().comm().allReduceInPlace(countNotInMask, eckit::mpi::sum());

      // Print stats about weights
      oops::Log::test() << "Mask " << maskParams.suffix.value() << " at level " << jlevel
        << ": inside/outside = " << countInMask << "/" << countNotInMask << std::endl;
    }

    // Add weight field
    wgtFset_->add(wgtField);
  }

  // Set name
  wgtFset_->name() = "weight";

  if (params_.filterParams.value()) {
    oops::Log::info() << "Info     : Filtering weight" << std::endl;

    // Create configuration without tests for filter
    const ErrorCovarianceParametersBase defaultParamsBase;

    // Create weight variables
    oops::Variables wgtVars;
    eckit::LocalConfiguration varConf;
    varConf.set("levels", maskField.levels());
    for (const auto & maskParams : params_.masksParams.value()) {
      wgtVars.push_back(oops::Variable(maskParams.suffix.value(), varConf));
    }

    // Initialize filter outer block chain
    const OuterBlockChain filter(outerGeometryData(),
                                 wgtVars,
                                 *fset4dXb_,
                                 *fset4dFg_,
                                 defaultParamsBase.toConfiguration(),
                                 *params_.filterParams.value());

    // Filter weight
    oops::FieldSet4D fset4d(*wgtFset_);
    filter.applyOuterBlocks(fset4d);
    wgtFset_->deepCopy(fset4d[0]);
  }

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  if (params_.inputWgtAtlasFileConf.value()) {
    // Create weight variables
    ASSERT(params_.maskLevels.value());
    oops::Variables wgtVars;
    eckit::LocalConfiguration varConf;
    varConf.set("levels", *params_.maskLevels.value());
    for (const auto & maskParams : params_.masksParams.value()) {
      wgtVars.push_back(oops::Variable(maskParams.suffix.value(), varConf));
    }

    // Read ATLAS weight file
    wgtFset_->read(innerGeometryData_.functionSpace(),
                   wgtVars,
                   *params_.inputWgtAtlasFileConf.value());
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>>
  GeographicalMask::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs;
  if (params_.outputWgtModelFileConf.value()) {
    outputs.push_back(std::make_pair(*params_.outputWgtModelFileConf.value(), *wgtFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void GeographicalMask::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Write ATLAS weight file
  if (params_.outputWgtAtlasFileConf.value()) {
    // Create weight variables
    oops::Variables wgtVars;
    eckit::LocalConfiguration varConf;
    for (const auto & field : *wgtFset_) {
      varConf.set("levels", field.levels());
      wgtVars.push_back(oops::Variable(field.name(), varConf));
    }

    // Print FieldSet norm
    oops::Log::test() << "Norm of output parameter weight: "
                      << wgtFset_->norm(wgtVars)
                      << std::endl;

    // Write file
    wgtFset_->write(*params_.outputWgtAtlasFileConf.value());
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}


// -----------------------------------------------------------------------------

void GeographicalMask::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
