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

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GeographicalMask> makerGeographicalMask_("GeographicalMask");

// -----------------------------------------------------------------------------

GeographicalMask::GeographicalMask(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerGeometryData_(outerGeometryData),
    activeOuterVars_(params.activeVars.value().get_value_or(outerVars)),
    params_(params.calibration.value() != boost::none ? *params.calibration.value()
      : *params.read.value()) {
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
  if (params_.inputMaskFileConf.value() && params_.inputWgtFileConf.value()) {
    throw eckit::UserError("cannot read both mask and weight at the same time", Here());
  }

  if (params.calibration.value()) {
    // Save configuration to read mask
    ASSERT(params_.inputMaskFileConf.value());
    readMaskConf_ = *params_.inputMaskFileConf.value();

    // Create mask fieldset
    maskFset_.reset(new oops::FieldSet3D(xb.validTime(), innerGeometryData_.comm()));
  }

  if (params.read.value()) {
    // Save configuration to read weight
    ASSERT(params_.inputWgtFileConf.value());
    readWgtConf_ = *params_.inputWgtFileConf.value();
  }

  // Create fieldsets
  wgtFset_.reset(new oops::FieldSet3D(xb.validTime(), innerGeometryData_.comm()));

  if (params.calibration.value() && params_.outputWgtFileConf.value()) {
    // Save configuration to write weight
    writeWgtConf_ = *params_.outputWgtFileConf.value();
  }

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
  if (!readMaskConf_.empty()) {
    inputs.push_back(std::make_pair("mask", readMaskConf_));
  }
  if (!readWgtConf_.empty()) {
    inputs.push_back(std::make_pair("weight", readWgtConf_));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void GeographicalMask::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  if (!readMaskConf_.empty()) {
    // Get mask from input file
    ASSERT(fsetVec.size() == 1);
    maskFset_->deepCopy(fsetVec[0]);
  }
  if (!readWgtConf_.empty()) {
    // Get weight from input file
    ASSERT(fsetVec.size() == 1);
    wgtFset_->deepCopy(fsetVec[0]);
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::directCalibration(const oops::FieldSets &) {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;

  // Create pointers
  ASSERT(maskFset_);
  ASSERT(wgtFset_);

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

    // TODO(Benjamin): smooth out weight if requested in yaml

    // Add wgt field
    wgtFset_->add(wgtField);
  }

  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeographicalMask::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  ASSERT(wgtFset_);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>>
  GeographicalMask::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs;

  if (!writeWgtConf_.empty()) {
    // Add pair
    outputs.push_back(std::make_pair(writeWgtConf_, *wgtFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void GeographicalMask::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
