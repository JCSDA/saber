/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberBlockChain.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberBlockChain::SaberBlockChain(const oops::Variables & incVars,
                                 const atlas::FieldSet & fset) :
  centralBlock_(), outerBlocks_(), centralFieldSet_(), incVars_(incVars), scalarWeightSqrt_(1.0),
  fileWeightSqrt_() {
  // Initialize central FieldSet
  centralFieldSet_ = util::copyFieldSet(fset);
  util::zeroFieldSet(centralFieldSet_);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::centralBlockInit(const oops::GeometryData & geometryData,
                                       const oops::Variables & outerVars,
                                       const eckit::Configuration & covarConf,
                                       const SaberBlockParametersBase & params,
                                       const atlas::FieldSet & xb,
                                       const atlas::FieldSet & fg,
                                       const util::DateTime & validTime,
                                       const size_t & timeRank) {
  centralBlock_ = SaberCentralBlockFactory::create(geometryData,
                                                   outerVars,
                                                   covarConf,
                                                   params,
                                                   xb,
                                                   fg,
                                                   validTime,
                                                   timeRank);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::setWeight(const double & scalarWeight) {
  ASSERT(scalarWeight > 0.0);
  scalarWeightSqrt_ = std::sqrt(scalarWeight);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::setWeight(const atlas::FieldSet & fileWeight) {
  fileWeightSqrt_ = util::copyFieldSet(fileWeight);
  util::sqrtFieldSet(fileWeightSqrt_);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::applyWeight(atlas::FieldSet & fset) const {
  // Weight square-root multiplication
  if (fileWeightSqrt_.empty()) {
    // Scalar weight
    util::multiplyFieldSet(fset, scalarWeightSqrt_);
  } else {
    // File-based weight
    util::multiplyFieldSets(fset, fileWeightSqrt_);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::applyOuterBlocks(atlas::FieldSet & fset) const {
  // Outer blocks forward multiplication
  for (auto it = outerBlocks_.rbegin(); it != outerBlocks_.rend(); ++it) {
    it->get()->multiply(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::applyOuterBlocksAD(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
    it->get()->multiplyAD(fset);
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::randomize(atlas::FieldSet & fset) const {
  // Initialize dx.fieldSet() with all the required fields, set to zero
  fset = util::copyFieldSet(centralFieldSet_);

  // Central block randomization
  centralBlock_->randomize(fset);

  // Outer blocks forward multiplication
  applyOuterBlocks(fset);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::multiply(atlas::FieldSet & fset) const {
  // Outer blocks adjoint multiplication
  applyOuterBlocksAD(fset);

  // Central block multiplication
  centralBlock_->multiply(fset);

  // Outer blocks forward multiplication
  applyOuterBlocks(fset);
}

// -----------------------------------------------------------------------------

void SaberBlockChain::leftInverseMultiply(atlas::FieldSet & fset) const {
  // Outer blocks left inverse multiplication
  for (auto it = outerBlocks_.begin(); it != outerBlocks_.end(); ++it) {
    if (it->get()->skipInverse()) {
      oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                      << it->get()->blockName() << std::endl;
    } else {
      it->get()->leftInverseMultiply(fset);
    }
  }
}

// -----------------------------------------------------------------------------

void SaberBlockChain::leftInverseMultiplyExceptLast(atlas::FieldSet & fset) const {
  // Outer blocks left inverse multiplication
  for (auto it = outerBlocks_.begin(); it != std::prev(outerBlocks_.end()); ++it) {
    if (it->get()->skipInverse()) {
      oops::Log::info() << "Warning: left inverse multiplication skipped for block "
                      << it->get()->blockName() << std::endl;
    } else {
      it->get()->leftInverseMultiply(fset);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
