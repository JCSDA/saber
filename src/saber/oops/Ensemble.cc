/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/Ensemble.h"

#include <vector>

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<Ensemble> makerEnsemble_("Ensemble");

// -----------------------------------------------------------------------------

Ensemble::Ensemble(const oops::GeometryData & geometryData,
                   const std::vector<size_t> & activeVariableSizes,
                   const oops::Variables & activeVars,
                   const Parameters_ & params,
                   const atlas::FieldSet & xb,
                   const atlas::FieldSet & fg,
                   const std::vector<atlas::FieldSet> & fsetVec)
{
  oops::Log::trace() << classname() << "::Ensemble starting" << std::endl;

  // Initialize ensemble
  for (const auto & fset : fsetVec) {
    ensemble_.push_back(fset);
  }

  // Initialize localization
  SaberCentralBlockParametersWrapper locParams;
  locParams.validateAndDeserialize(params.localization);
  const SaberBlockParametersBase & saberCentralBlockParams =
    locParams.saberCentralBlockParameters;

  oops::Log::info() << "Info     : Creating localization block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;
  loc_.reset(SaberCentralBlockFactory::create(geometryData, activeVariableSizes, activeVars,
    saberCentralBlockParams, xb, fg, fsetVec));

  oops::Log::trace() << classname() << "::Ensemble done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Initialization
  util::zeroFieldSet(fset);

  // Loop over ensemble members
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

  // Loop over ensemble members
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

  // Normalize result
  const double rk = 1.0/static_cast<double>(ensemble_.size()-1);
  util::multiplyFieldSet(fset, rk);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Ensemble::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
