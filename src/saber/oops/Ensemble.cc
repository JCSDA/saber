/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/Ensemble.h"

#include <utility>
#include <vector>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "saber/oops/SaberBlockChain.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<Ensemble> makerEnsemble_("Ensemble");

// -----------------------------------------------------------------------------

Ensemble::Ensemble(const oops::GeometryData & geometryData,
                   const std::vector<size_t> & activeVariableSizes,
                   const oops::Variables & activeVars,
                   const eckit::Configuration & covarConf,
                   const Parameters_ & params,
                   const atlas::FieldSet & xb,
                   const atlas::FieldSet & fg,
                   const util::DateTime & validTimeOfXbFg,
                   const size_t & timeRank) :
  SaberCentralBlockBase(params),
  timeRank_(timeRank), vars_(activeVars), comm_(geometryData.comm()) {
  oops::Log::trace() << classname() << "::Ensemble starting" << std::endl;
  oops::Log::trace() << classname() << "::Ensemble done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Initialization
  util::zeroFieldSet(fset);

  if (loc_) {
    // Apply localization
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Temporary copy for this ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(fset);

      // Randomize SABER central block
      loc_->randomize(fsetMem);

      // Schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  } else {
    // No localization
    util::NormalDistribution<double> normalDist(ensemble_.size(), 0.0, 1.0, seed_);
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Copy ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(ensemble_[ie]);

      // Apply weight
      util::multiplyFieldSet(fsetMem, normalDist[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  }

  // Normalize result
  const double rk = 1.0/sqrt(static_cast<double>(ensemble_.size()-1));
  util::multiplyFieldSet(fset, rk);

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Initialization
  atlas::FieldSet fsetInit = util::copyFieldSet(fset);
  util::zeroFieldSet(fset);

  if (loc_) {
    // Apply localization
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Temporary copy for this ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(fsetInit);

      // First schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Apply localization
      loc_->multiply(fsetMem);

      // Second schur product
      util::multiplyFieldSets(fsetMem, ensemble_[ie]);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  } else {
    // No localization
    for (unsigned int ie = 0; ie < ensemble_.size(); ++ie) {
      // Compute weight
      const double wgt = util::dotProductFieldSets(fsetInit,
                                                   ensemble_[ie],
                                                   vars_.variables(),
                                                   comm_,
                                                   false);

      // Copy ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(ensemble_[ie]);

      // Apply weight
      util::multiplyFieldSet(fsetMem, wgt);

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }
  }

  // Normalize result
  const double rk = 1.0/static_cast<double>(ensemble_.size()-1);
  util::multiplyFieldSet(fset, rk);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::directCalibration(const std::vector<atlas::FieldSet> & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;
  // Initialize ensemble
  for (const auto & fset : fsetEns) {
    ensemble_.push_back(fset);
  }
  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::setLocalization(std::unique_ptr<SaberBlockChain> locBlockChain) {
  // Pass central block
  loc_ = std::move(locBlockChain);
}

// -----------------------------------------------------------------------------

void Ensemble::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
