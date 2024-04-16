/*
 * (C) Copyright 2024 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <omp.h>

#include <algorithm>
#include <string>

#include "saber/generic/FakeLevels.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<FakeLevels> makerFakeLevels_("FakeLevels");

// -----------------------------------------------------------------------------

FakeLevels::FakeLevels(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    validTime_(xb.validTime()),
    gdata_(outerGeometryData),
    comm_(gdata_.comm()),
    outerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars_)),
    suffix_("_fakeLevels"),
    params_(params.calibration.value() != boost::none ? *params.calibration.value()
      : *params.read.value()),
    fieldsMetaData_(params.fieldsMetaData.value()) {
  oops::Log::trace() << classname() << "::FakeLevels starting" << std::endl;

  // Set number of fake levels
  if (params_.fakeLevels.value() != boost::none) {
    nz_ = params_.fakeLevels.value()->size();
    if (params_.nz.value() != boost::none) {
      ASSERT(nz_ = *params_.nz.value());
    }
  } else if (params_.nz.value() != boost::none) {
    nz_ = *params_.nz.value();
  } else {
    throw eckit::UserError("number of fake levels is missing", Here());
  }
  ASSERT(nz_ > 1);

  // Create inner variables
  for (const std::string & varName : outerVars.variables()) {
    if (activeVars_.has(varName)) {
      const std::string newVarName = varName + suffix_;
      innerVars_.push_back(newVarName);
      innerVars_.addMetaData(newVarName, "levels", nz_);
    } else {
      innerVars_.push_back(varName);
      innerVars_.addMetaData(varName, "levels", innerVars_.getLevels(varName));
    }
  }

  oops::Log::trace() << classname() << "::FakeLevels done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet outerFset;

  for (const auto & outerVar : activeVars_.variables()) {
    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Get inner field
    const auto innerView = atlas::array::make_view<double, 2>(fset[innerVar]);

    // Create outer field
    atlas::Field outerField = gdata_.functionSpace().createField<double>(
      atlas::option::name(outerVar) | atlas::option::levels(1) | atlas::option::halo(1));
    auto outerView = atlas::array::make_view<double, 2>(outerField);

    // Get weight
    const auto weightView = atlas::array::make_view<double, 2>((*weight_)[innerVar]);

    // Reduce to a single level
    outerView.assign(0.0);
    for (int jnode = 0; jnode < outerField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          outerView(jnode, 0) += innerView(jnode, k)*weightView(jnode, k);
        }
      }
    }

    // Add field
    outerFset.add(outerField);
  }

  // Copy outer fieldset
  fset.fieldSet() = outerFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Create fieldset
  atlas::FieldSet innerFset;

  for (const auto & outerVar : activeVars_.variables()) {
    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Get outer field
    const auto outerView = atlas::array::make_view<double, 2>(fset[outerVar]);
    ASSERT(fset[outerVar].shape(1) == 1);

    // Create inner field
    atlas::Field innerField = gdata_.functionSpace().createField<double>(
      atlas::option::name(innerVar) | atlas::option::levels(nz_) | atlas::option::halo(1));
    auto innerView = atlas::array::make_view<double, 2>(innerField);

    // Get weight
    const auto weightView = atlas::array::make_view<double, 2>((*weight_)[innerVar]);

    // Extend to multiple levels
    innerView.assign(0.0);
    for (int jnode = 0; jnode < innerField.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t k = 0; k < nz_; ++k) {
          innerView(jnode, k) = outerView(jnode, 0)*weightView(jnode, k);
        }
      }
    }

    // Add field
    innerFset.add(innerField);
  }

  // Copy inner fieldset
  fset.fieldSet() = innerFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> FakeLevels::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (params_.inputModelFilesConf.value() != boost::none) {
    for (const auto & conf : *params_.inputModelFilesConf.value()) {
      // File name
      const std::string fileName = conf.getString("parameter");

      // Get file configuration
      eckit::LocalConfiguration fileConf = getFileConf(comm_, conf);

      // Add pair
      inputs.push_back(std::make_pair(fileName, fileConf));
    }
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void FakeLevels::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  // Get rv from input files if present
  for (size_t ji = 0; ji < fsetVec.size(); ++ji) {
    if (fsetVec[ji].name() == "rv") {
      rv_.reset(new oops::FieldSet3D(validTime_, comm_));
      for (const auto & outerVar : activeVars_.variables()) {
        // Save field
        rv_->add(fsetVec[ji][outerVar]);
      }
    } else if (fsetVec[ji].name() == "weight") {
      weight_.reset(new oops::FieldSet3D(validTime_, comm_));
      for (const auto & outerVar : activeVars_.variables()) {
        // Get inner variable name
        const std::string innerVar = outerVar + suffix_;

        // Save field
        weight_->add(fsetVec[ji][innerVar]);
      }
    } else {
      throw eckit::UserError("wrong input parameter", Here());
    }
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::directCalibration(const oops::FieldSets &) {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Get vertical support from yaml if needed
  if (!rv_) {
    ASSERT(params_.rvFromYaml.value() != boost::none);
    rv_.reset(new oops::FieldSet3D(validTime_, comm_));
    for (const auto & outerVar : activeVars_.variables()) {
      atlas::Field rvField = gdata_.functionSpace().createField<double>(
        atlas::option::name(outerVar) | atlas::option::levels(1));
      auto rvView = atlas::array::make_view<double, 2>(rvField);
      for (int jnode = 0; jnode < rvField.shape(0); ++jnode) {
        if (ghostView(jnode) == 0) {
          rvView(jnode, 0) = *params_.rvFromYaml.value();
        }
      }
      rv_->add(rvField);
    }
  }

  // Define fake levels
  std::vector<double> fakeLevels;
  if (params_.fakeLevels.value() != boost::none) {
    fakeLevels = *params_.fakeLevels.value();
  } else {
    ASSERT(params_.lowestFakeLevel.value() != boost::none);
    ASSERT(params_.highestFakeLevel.value() != boost::none);
    ASSERT(*params_.lowestFakeLevel.value() < *params_.highestFakeLevel.value());
    const double delta = (*params_.highestFakeLevel.value()-*params_.lowestFakeLevel.value())
      /(nz_-1);
    for (size_t k = 0; k < nz_; ++k) {
      fakeLevels.push_back(*params_.lowestFakeLevel.value()+static_cast<double>(k)*delta);
    }
  }

  // Prepare weight
  ASSERT(!weight_);
  weight_.reset(new oops::FieldSet3D(validTime_, comm_));
  for (const auto & outerVar : activeVars_.variables()) {
    // Check number of levels of outer variables
    ASSERT(outerVars_.getLevels(outerVar) == 1);

    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Get vertical coordinate
    const std::string key = outerVar + ".vert_coord";
    const std::string vertCoordName = fieldsMetaData_.getString(key, "vert_coord");
    const atlas::Field vertCoordField = gdata_.fieldSet()[vertCoordName];
    const auto vertCoordView = atlas::array::make_view<double, 2>(vertCoordField);

    // Get vertical support
    const auto rvView = atlas::array::make_view<double, 2>((*rv_)[outerVar]);

    // Create weight field
    atlas::Field field = gdata_.functionSpace().createField<double>(
      atlas::option::name(innerVar) | atlas::option::levels(nz_) | atlas::option::halo(1));

    // Set to zero
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(0.0);

    // Compute weight
    std::vector<double> wgt(nz_);
    for (int jnode = 0; jnode < field.shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Check fake levels extrema
        ASSERT(vertCoordView(jnode, 0) >= fakeLevels[0]);
        ASSERT(vertCoordView(jnode, 0) <= fakeLevels[nz_-1]);

        // Compute raw weight
        double wgtSum = 0.0;
        for (size_t k = 0; k < nz_; ++k) {
          const double normDist = std::abs(vertCoordView(jnode, 0)-fakeLevels[k])/rvView(jnode, 0);
          if ((fakeLevels[k] < vertCoordView(jnode, 0)) && (normDist > 1.0e-12)) {
            // Under ground
            wgt[k] = 0.0;
          } else {
            // Above ground
            wgt[k] = oops::gc99(normDist);
            wgt[k] = std::max(0.0, wgt[k]);
          }
          wgtSum += wgt[k];
        }

        // Normalize weight
        if (wgtSum > 0) {
          for (size_t k = 0; k < nz_; ++k) {
            wgt[k] /= wgtSum;
          }
        } else {
          throw eckit::UserError("fake levels are too far apart", Here());
        }

        // Weight square-roots
        for (size_t k = 0; k < nz_; ++k) {
          ASSERT(wgt[k] >= 0.0);
          wgt[k] = std::sqrt(wgt[k]);
        }

        // Copy weight
        for (size_t k = 0; k < nz_; ++k) {
          view(jnode, k) = wgt[k];
        }
      }
    }

    // Add field
    weight_->add(field);
  }

  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void FakeLevels::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  ASSERT(!rv_);
  ASSERT(weight_);
  for (const auto & outerVar : activeVars_.variables()) {
    // Get inner variable name
    const std::string innerVar = outerVar + suffix_;

    // Check number of levels
    ASSERT((*weight_)[innerVar].levels() == static_cast<int>(nz_));
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> FakeLevels::fieldsToWrite()
  const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> pairs;

  // Get vector of configurations
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = params_.outputModelFilesConf.value().get_value_or({});

  for (const auto & conf : outputModelFilesConf) {
    // Get file configuration
    eckit::LocalConfiguration file = getFileConf(comm_, conf);

    // Get parameter
    const std::string param = conf.getString("parameter");

    // Add pair
    if (param == "rv") {
      rv_->name() = "rv";
      pairs.push_back(std::make_pair(file, *rv_));
    } else if (param == "weight") {
      weight_->name() = "weight";
      pairs.push_back(std::make_pair(file, *weight_));
    } else {
      throw eckit::UserError("wrong output parameter", Here());
    }
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void FakeLevels::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration FakeLevels::getFileConf(const eckit::mpi::Comm & comm,
                                                  const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::getFileConf starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Get IO configuration
  eckit::LocalConfiguration file = conf.getSubConfiguration("file");

  // Replace patterns
  util::seekAndReplace(file, "_MPI_", mpi);
  util::seekAndReplace(file, "_OMP_", omp);

  oops::Log::trace() << classname() << "::getFileConf done" << std::endl;
  return file;
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
