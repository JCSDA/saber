/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/StdDev.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<StdDev> makerStdDev_("StdDev");

// -----------------------------------------------------------------------------

StdDev::StdDev(const oops::GeometryData & outerGeometryData,
               const oops::Variables & outerVars,
               const eckit::Configuration & covarConfig,
               const Parameters_ & params,
               const oops::FieldSet3D & xb,
               const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    params_(params),
    readFromAtlas_(false),
    readFromModel_(false),
    scaleFactor_(params_.scaleFactorParam.value()),
    stdDevFset_(),
    writeToAtlas_(false),
    writeToModel_(false)
{
  oops::Log::trace() << classname() << "::StdDev starting" << std::endl;

  // Prepare read parameters
  const auto & readParams = params_.readParams.value();
  if (readParams != boost::none) {
    const auto & atlasFileConf = readParams->atlasFileConf.value();
    readFromAtlas_ = (atlasFileConf != boost::none);
    const auto & modelFileConf = readParams->modelFileConf.value();
    readFromModel_ = (modelFileConf != boost::none);
    if (readFromAtlas_ && readFromModel_) {
      throw eckit::UserError("either ATLAS or model stddev input file should be present, not both",
        Here());
    }
    if (readFromAtlas_) {
      readConf_ = *atlasFileConf;
    }
    if (readFromModel_) {
      readConf_ = *modelFileConf;
    }
    stdDevFset_.reset(new oops::FieldSet3D(xb.validTime(), innerGeometryData_.comm()));
  }

  // Prepare write parameters
  const auto & calibrationParams = params_.calibrationParams.value();
  if (calibrationParams != boost::none) {
    const auto & atlasFileConf = calibrationParams->atlasFileConf.value();
    writeToAtlas_ = (atlasFileConf != boost::none);
    const auto & modelFileConf = calibrationParams->modelFileConf.value();
    writeToModel_ = (modelFileConf != boost::none);
    if (writeToAtlas_ && writeToModel_) {
      throw eckit::UserError("either ATLAS or model stddev output file should be present, not both",
        Here());
    }
    if (writeToAtlas_) {
      writeConf_ = *atlasFileConf;
    }
    if (writeToModel_) {
      writeConf_ = *modelFileConf;
    }
  }

  const eckit::LocalConfiguration & scaling = params_.scaling.value();
  std::vector<eckit::LocalConfiguration> scales = scaling.getSubConfigurations();
  for (const eckit::LocalConfiguration & scale : scales) {
    const std::string var = scale.getString("variable");
    scaling_[var] = scale.getDouble("stddev");
  }

  oops::Log::trace() << classname() << "::StdDev done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  if (stdDevFset_) {
    fset *= *stdDevFset_;
  } else {
    for (auto & field : fset) {
      const std::string var = field.name();
      if (scaling_.find(var) != scaling_.end()) {
        const double fact = scaling_.at(var);
        auto view = atlas::array::make_view<double, 2>(field);
        for (int jnode = 0; jnode < field.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) *= fact;
          }
        }
      }
    }
  }
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  if (stdDevFset_) {
    fset *= *stdDevFset_;
  } else {
    for (auto & field : fset) {
      const std::string var = field.name();
      if (scaling_.find(var) != scaling_.end()) {
        const double fact = scaling_.at(var);
        auto view = atlas::array::make_view<double, 2>(field);
        for (int jnode = 0; jnode < field.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) *= fact;
          }
        }
      }
    }
  }
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  if (stdDevFset_) {
    fset /= *stdDevFset_;
  } else {
    for (auto & field : fset) {
      const std::string var = field.name();
      if (scaling_.find(var) != scaling_.end()) {
        const double fact = 1.0 / scaling_.at(var);
        auto view = atlas::array::make_view<double, 2>(field);
        for (int jnode = 0; jnode < field.shape(0); ++jnode) {
          for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
            view(jnode, jlevel) *= fact;
          }
        }
      }
    }
  }
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> StdDev::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (readFromModel_) {
    inputs.push_back(std::make_pair("StdDev", readConf_));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void StdDev::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  if (readFromModel_) {
    ASSERT(fsetVec.size() == 1);
    stdDevFset_->deepCopy(fsetVec[0]);

    // apply scale factor (default = 1.0)
    *stdDevFset_ *= this->scaleFactor_;
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Read ATLAS stddev file
  if (readFromAtlas_) {
    // Read file
    stdDevFset_->read(innerGeometryData_.functionSpace(),
                      innerVars_,
                      readConf_);

    // apply scale factor (default = 1.0)
    *stdDevFset_ *= this->scaleFactor_;

    // Set name
    stdDevFset_->name() = "StdDev";

    // Print FieldSet norm
    oops::Log::test() << "Norm of input parameter StdDev: "
                      << stdDevFset_->norm(innerVars_)
                      << std::endl;
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void StdDev::directCalibration(const oops::FieldSets & fsetEns) {
  // Initialize
  oops::FieldSet3D mean(this->validTime(), innerGeometryData_.comm());
  oops::FieldSet3D var(this->validTime(), innerGeometryData_.comm());
  for (const auto & innerVar : innerVars_) {
    mean.fieldSet().add(innerGeometryData_.functionSpace().createField<double>(
             atlas::option::name(innerVar.name()) |
             atlas::option::levels(innerVar.getLevels())));
    var.fieldSet().add(innerGeometryData_.functionSpace().createField<double>(
            atlas::option::name(innerVar.name()) |
            atlas::option::levels(innerVar.getLevels())));
  }
  mean.zero();
  var.zero();

  // Compute mean
  for (size_t jj = 0; jj < fsetEns.ens_size(); ++jj) {
    mean += fsetEns[jj];
  }
  mean *= 1.0/static_cast<double>(fsetEns.ens_size());

  // Compute variance
  for (size_t jj = 0; jj < fsetEns.ens_size(); ++jj) {
    // Compute perturbation
    oops::FieldSet3D fset_pert(fsetEns[jj]);
    fset_pert -= mean;

    // Compute squared pertrubation
    fset_pert *= fset_pert;

    // Update sum
    var += fset_pert;
  }
  var *= 1.0/static_cast<double>(fsetEns.ens_size()-1);

  // Get standard deviation
  stdDevFset_.reset(new oops::FieldSet3D(var));
  stdDevFset_->sqrt();

  // Set name
  stdDevFset_->name() = "StdDev";
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationInit() {
  // Initialize iterative counters with zeroes
  iterativeN_ = 0;
  iterativeMean_.reset(new oops::FieldSet3D(this->validTime(), innerGeometryData_.comm()));
  iterativeVar_.reset(new oops::FieldSet3D(this->validTime(), innerGeometryData_.comm()));
  for (const auto & innerVar : innerVars_) {
    iterativeMean_->fieldSet().add(innerGeometryData_.functionSpace().createField<double>(
                       atlas::option::name(innerVar.name()) |
                       atlas::option::levels(innerVar.getLevels())));
    iterativeVar_->fieldSet().add(innerGeometryData_.functionSpace().createField<double>(
                      atlas::option::name(innerVar.name()) |
                      atlas::option::levels(innerVar.getLevels())));
  }
  iterativeMean_->zero();
  iterativeVar_->zero();
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  // Increment ensemble index
  // ie = ie + 1
  iterativeN_++;

  // Remove mean
  // pert = state - mean
  oops::FieldSet3D fset_pert(fset);
  fset_pert -= *iterativeMean_;

  // Update variance
  // var = var + (ie-1)/ie * pert^2
  oops::FieldSet3D fset_pert_squared(fset_pert);
  fset_pert_squared *= fset_pert;
  fset_pert_squared *= static_cast<double>(iterativeN_-1)/static_cast<double>(iterativeN_);
  *iterativeVar_ += fset_pert_squared;

  // Update mean
  // mean = mean + 1 / ie * pert
  fset_pert *= 1.0/static_cast<double>(iterativeN_);
  *iterativeMean_ += fset_pert;
}

// -----------------------------------------------------------------------------

void StdDev::iterativeCalibrationFinal() {
  // Normalize variance
  // var = 1 / (N-1) * var
  *iterativeVar_ *= 1.0/static_cast<double>(iterativeN_-1);

  // Get standard-deviation
  stdDevFset_.reset(new oops::FieldSet3D(*iterativeVar_));
  stdDevFset_->sqrt();

  // Set name
  stdDevFset_->name() = "StdDev";

  // Cleaning
  iterativeN_ = 0;
  iterativeMean_.reset();
  iterativeVar_.reset();
}

// -----------------------------------------------------------------------------

void StdDev::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Write ATLAS stddev file
  if (writeToAtlas_) {
    // Print FieldSet norm
    oops::Log::test() << "Norm of output parameter StdDev: "
                      << stdDevFset_->norm(innerVars_)
                      << std::endl;

    // Write file
    stdDevFset_->write(writeConf_);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> StdDev::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs;

  if (writeToModel_) {
    // Add pair
    outputs.push_back(std::make_pair(writeConf_, *stdDevFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return outputs;
}

// -----------------------------------------------------------------------------

void StdDev::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
