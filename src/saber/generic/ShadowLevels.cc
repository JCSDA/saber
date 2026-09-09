/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <omp.h>

#include <algorithm>
#include <limits>
#include <string>

#include "saber/generic/ShadowLevels.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<ShadowLevels> makerShadowLevels_("ShadowLevels");

// -----------------------------------------------------------------------------

ShadowLevels::ShadowLevels(const oops::GeometryData & outerGeometryData,
                           const oops::Variables & outerVars,
                           const eckit::Configuration & covarConf,
                           const Parameters_ & params,
                           const oops::FieldSet3D & xb,
                           const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    gdata_(outerGeometryData),
    comm_(gdata_.comm()),
    activeOuterVars_(params.activeVars.value().get_value_or(outerVars)),
    params_(params.calibration.value() != boost::none ? *params.calibration.value()
      : *params.read.value()),
    fieldsMetaData_(params.fieldsMetaData.value()) {
  oops::Log::trace() << classname() << "::ShadowLevels starting" << std::endl;

  // Check that all active variables are 2D
  for (const auto & var : activeOuterVars_) {
    ASSERT(var.getLevels() == 1);
  }

  // Copy groups
  std::vector<std::string> outerVarsCheck;
  for (const auto & groupParams : params_.groups.value()) {
    // Define group properties
    Group group;
    group.suffix_ = groupParams.suffix.value();
    if (groupParams.name.value()) {
      group.name_ = *groupParams.name.value();
    } else {
      ASSERT(groupParams.variables.value().size() == 1);
      group.name_ = groupParams.variables.value()[0] + group.suffix_;
    }
    group.variables_ = groupParams.variables.value();
    group.nz_ = groupParams.nz.value();
    group.rv_ = groupParams.rv.value() ? *groupParams.rv.value() : -1.0;
    group.cmpWgtSqrt_ = std::sqrt(groupParams.cmpWgt.value());

    // Check variables consistency
    for (const auto & varName : group.variables_) {
      outerVarsCheck.push_back(varName);
    }

    // Add group
    groups_.push_back(group);
  }

  // Check that active variables are all present in groups
  for (const auto & var : activeOuterVars_) {
    ASSERT(std::find(outerVarsCheck.begin(), outerVarsCheck.end(), var.name())
      != outerVarsCheck.end());
  }

  // Add active inner variables
  for (const auto & group : groups_) {
    for (const auto & varName : group.variables_) {
      const std::string innerVarName = varName + group.suffix_;
      ASSERT(!innerVars_.has(innerVarName));
      eckit::LocalConfiguration conf;
      conf.set("levels", group.nz_);
      innerVars_.push_back(oops::Variable(innerVarName, conf));
    }
  }

  // Add inactive inner variables
  for (const auto & outerVar : outerVars) {
    if (!activeOuterVars_.has(outerVar)) {
      innerVars_.push_back(outerVars[outerVar.name()]);
    }
  }

  if (params.read.value()) {
    // Save configuration to read weight
    ASSERT(params_.inputWgtFileConf.value());
    readConf_ = *params_.inputWgtFileConf.value();

    // Create mask fieldset
    wgtFset_.reset(new oops::FieldSet3D(xb.validTime(), gdata_.comm()));
  }

  if (params.calibration.value() && params_.outputWgtFileConf.value()) {
    // Save configuration to write weight
    writeConf_ = *params_.outputWgtFileConf.value();
  }

  oops::Log::trace() << classname() << "::ShadowLevels done" << std::endl;
}

// -----------------------------------------------------------------------------

void ShadowLevels::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet outerFset;

  for (const auto & group : groups_) {
    // Get weight view
    const auto weightView = atlas::array::make_view<double, 2>((*wgtFset_)[group.name_]);

    for (const auto & varName : group.variables_) {
      // Get inner variable name
      const std::string innerVarName = varName + group.suffix_;

      // Get inner field
      const auto innerView = atlas::array::make_view<double, 2>(fset[innerVarName]);

      // Outer field
      if (!outerFset.has(varName)) {
        // Create outer field
        atlas::Field outerField = gdata_.functionSpace().createField<double>(
          atlas::option::name(varName) | atlas::option::levels(1) | atlas::option::halo(1));

        // Get outer view
        auto outerView = atlas::array::make_view<double, 2>(outerField);

        // Initialize outer view
        outerView.assign(0.0);

        // Add outer field
        outerFset.add(outerField);
      }

      // Get outer field
      atlas::Field outerField = outerFset[varName];

      // Get outer field view
      auto outerView = atlas::array::make_view<double, 2>(outerField);

      // Reduce to a single level
      for (int jnode = 0; jnode < outerField.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          for (size_t k = 0; k < group.nz_; ++k) {
            outerView(jnode, 0) += innerView(jnode, k)*weightView(jnode, k);
          }
        }
      }
    }
  }

  // Add inactive variables
  for (const auto & outerVar : outerVars()) {
    if (!activeOuterVars_.has(outerVar.name())) {
      outerFset.add(fset[outerVar.name()]);
    }
  }

  // Copy outer fieldset
  fset.fieldSet() = outerFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ShadowLevels::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet innerFset;

  for (const auto & group : groups_) {
    // Get weight
    const auto weightView = atlas::array::make_view<double, 2>((*wgtFset_)[group.name_]);

    for (const auto & varName : group.variables_) {
      // Get inner variable name
      const std::string innerVarName = varName + group.suffix_;

      // Get outer field
      const auto outerView = atlas::array::make_view<double, 2>(fset[varName]);
      ASSERT(fset[varName].shape(1) == 1);

      // Create inner field
      atlas::Field innerField = gdata_.functionSpace().createField<double>(
        atlas::option::name(innerVarName) | atlas::option::levels(group.nz_)
        | atlas::option::halo(1));
      auto innerView = atlas::array::make_view<double, 2>(innerField);

      // Extend to multiple levels
      innerView.assign(0.0);
      for (int jnode = 0; jnode < innerField.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          for (size_t k = 0; k < group.nz_; ++k) {
            innerView(jnode, k) = outerView(jnode, 0)*weightView(jnode, k);
          }
        }
      }

      // Add field
      innerFset.add(innerField);
    }
  }

  // Add inactive variables
  for (const auto & outerVar : outerVars()) {
    if (!activeOuterVars_.has(outerVar.name())) {
      innerFset.add(fset[outerVar.name()]);
    }
  }

  // Copy inner fieldset
  fset.fieldSet() = innerFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> ShadowLevels::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (!readConf_.empty()) {
    inputs.push_back(std::make_pair("weight", readConf_));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void ShadowLevels::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  if (!readConf_.empty()) {
    // Get weight from input file
    ASSERT(fsetVec.size() == 1);
    wgtFset_->deepCopy(fsetVec[0]);
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void ShadowLevels::directCalibration(const oops::FieldSets &) {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;

  // Prepare weight
  ASSERT(!wgtFset_);
  wgtFset_.reset(new oops::FieldSet3D(validTime_, comm_));

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  for (const auto & group : groups_) {
    // Print group properties
    oops::Log::test() << "Group: " << group.name_ << std::endl;
    oops::Log::test() << "- variables: " << group.variables_ << std::endl;
    oops::Log::test() << "- rv: " << group.rv_ << std::endl;

    // Check rv value
    ASSERT(group.rv_ > 0.0);

    // Apply GC99 autoconvolution factor
    const double rv = group.rv_*autoConvolFactor_;

    // Equivalent Gaussian length-scale
    const double L = rv*rv2L_;

    // Get vertical coordinate for this group
    const std::string key = group.variables_[0] + ".vert_coord";
    const std::string vertCoordName = fieldsMetaData_.getString(key, "vert_coord");
    const atlas::Field vertCoordField = gdata_.fieldSet()[vertCoordName];
    const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoordField);

    // Get min and max of the vertical coordinate
    double vcMin = std::numeric_limits<double>::max();
    double vcMax = std::numeric_limits<double>::min();
    for (size_t jnode = 0; jnode < vertCoordField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        vcMin = std::min(vcMin, vertCoordView(jnode, 0));
        vcMax = std::max(vcMax, vertCoordView(jnode, 0));
      }
    }
    comm_.allReduceInPlace(vcMin, eckit::mpi::min());
    comm_.allReduceInPlace(vcMax, eckit::mpi::max());
    oops::Log::test() << "- min/max vertical coordinate: " << vcMin << " / " << vcMax << std::endl;

    // Vertical coordinate span
    const double vcSpan = vcMax-vcMin;

    // SABER assumption: maximum shadow levels spacing dh = L, extension layer thickness de = 0

    // Minimum number of shadow levels
    const size_t nz = static_cast<size_t>(vcSpan/L)+2;
    oops::Log::test() << "- minimum number of shadow levels: " << nz << std::endl;
    oops::Log::test() << "- actual number of shadow levels: " << group.nz_ << std::endl;

    if (group.nz_ < nz) {
      // Not enough shadow levels
      const std::string message = "Not enough shadow levels, at least " + std::to_string(nz)
        + " levels required";
      throw eckit::UserError(message, Here());
    }

    // Set shadow levels spacing
    const double dh = vcSpan/static_cast<double>(group.nz_-1);
    oops::Log::test() << "- normalized thickness: " << dh/L << std::endl;

    // Define shadow levels
    std::vector<double> shadowLevels;
    for (size_t k = 0; k < group.nz_; ++k) {
      shadowLevels.push_back(vcMin+static_cast<double>(k)*dh);
    }
    ASSERT(oops::is_close_relative(shadowLevels[group.nz_-1], vcMax, 1.0e-12));

    // Create weight field
    atlas::Field field = gdata_.functionSpace().createField<double>(
      atlas::option::name(group.name_) | atlas::option::levels(group.nz_)
      | atlas::option::halo(1));

    // Set to zero
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);

    // Compute weight
    std::vector<double> wgt(group.nz_);
    std::vector<double> sqWgtAvg(group.nz_, 0.0);
    int numberOfCells = 0;
    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Compute raw weight
        double wgtSum = 0.0;
        for (size_t k = 0; k < group.nz_; ++k) {
          const double normDist = std::abs(vertCoordView(jnode, 0)-shadowLevels[k])/rv;
          wgt[k] = oops::gc99(normDist);
          wgt[k] = std::max(0.0, wgt[k]);
          wgtSum += wgt[k];
        }

        // Normalize weight
        if (wgtSum > 0) {
          for (size_t k = 0; k < group.nz_; ++k) {
            wgt[k] /= wgtSum;
          }
        } else {
          throw eckit::UserError("shadow levels are too far apart", Here());
        }

        // Weight square-roots
        for (size_t k = 0; k < group.nz_; ++k) {
          ASSERT(wgt[k] >= 0.0);
          wgt[k] = std::sqrt(wgt[k]);
        }

        // Apply component weight
        for (size_t k = 0; k < group.nz_; ++k) {
          view(jnode, k) = group.cmpWgtSqrt_*wgt[k];
        }

        // Horizontally averaged squared weight (for diagnostic purpose)
        for (size_t k = 0; k < group.nz_; ++k) {
          sqWgtAvg[k] += wgt[k]*wgt[k];
        }
        ++numberOfCells;
      }
    }

    // Add field
    wgtFset_->add(field);

    // Reduce averaged weight
    comm_.allReduceInPlace(sqWgtAvg.begin(), sqWgtAvg.end(), eckit::mpi::sum());
    comm_.allReduceInPlace(numberOfCells, eckit::mpi::sum());
    const std::ios_base::fmtflags testFlags = oops::Log::test().flags();
    const std::streamsize testPrecision = oops::Log::test().precision();
    oops::Log::test() << "- squared weight average: " << std::fixed << std::setprecision(1);
    for (const auto & lev : sqWgtAvg) {
      oops::Log::test() << 100.0*lev/static_cast<double>(numberOfCells) << "% ";
    }
    oops::Log::test().flags(testFlags);
    oops::Log::test().precision(testPrecision);
    oops::Log::test() << std::endl;
  }

  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void ShadowLevels::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  ASSERT(wgtFset_);
  for (auto & group : groups_) {
    // Check number of levels
    const size_t nzInFile = (*wgtFset_)[group.name_].levels();
    if (nzInFile != group.nz_) {
      const std::string message = "Number of levels in file [" + std::to_string(nzInFile) + "] is "
        + "different from the number number of shadow levels [" + std::to_string(group.nz_) + "]";
      throw eckit::UserError(message, Here());
    }
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> ShadowLevels::fieldsToWrite()
  const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs;

  if (!writeConf_.empty()) {
    // Add pair
    outputs.push_back(std::make_pair(writeConf_, *wgtFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void ShadowLevels::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
