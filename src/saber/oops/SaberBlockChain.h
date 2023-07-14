/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberBlockChain {
 public:
  SaberBlockChain();
  explicit SaberBlockChain(const oops::Variables &,
                           const atlas::FieldSet &);
  ~SaberBlockChain() {}

  // Central block initialization
  void centralBlockInit(const oops::GeometryData &,
                        const oops::Variables &,
                        const eckit::Configuration &,
                        const SaberBlockParametersBase &,
                        const atlas::FieldSet &,
                        const atlas::FieldSet &,
                        const util::DateTime &,
                        const size_t & timeRank = 0);

  // Hybrid weight initialization
  void setWeight(const double &);
  void setWeight(const atlas::FieldSet &);
  void applyWeight(atlas::FieldSet &) const;

  // Accessors
  const SaberCentralBlockBase & centralBlock() const {return *centralBlock_;}
  SaberCentralBlockBase & centralBlock() {return *centralBlock_;}
  const std::vector<std::unique_ptr<SaberOuterBlockBase>> & outerBlocks() const
    {return outerBlocks_;}
  std::vector<std::unique_ptr<SaberOuterBlockBase>> & outerBlocks() {return outerBlocks_;}
  const SaberOuterBlockBase & lastOuterBlock() const {return *outerBlocks_.back();}
  SaberOuterBlockBase & lastOuterBlock() {return *outerBlocks_.back();}
  const atlas::FieldSet & centralFieldSet() const {return centralFieldSet_;}
  atlas::FieldSet & centralFieldSet() {return centralFieldSet_;}
  const oops::Variables & incVars() const {return incVars_;}

  // Apply outer blocks (forward and adjoint)
  void applyOuterBlocks(atlas::FieldSet &) const;
  void applyOuterBlocksAD(atlas::FieldSet &) const;

  // Randomize and multiply
  void randomize(atlas::FieldSet &) const;
  void multiply(atlas::FieldSet &) const;

  // Calibration inverse multiply
  void leftInverseMultiply(atlas::FieldSet &) const;
  void leftInverseMultiplyExceptLast(atlas::FieldSet &) const;

 private:
  // Central block
  std::unique_ptr<SaberCentralBlockBase> centralBlock_;

  // Outer blocks
  std::vector<std::unique_ptr<SaberOuterBlockBase>> outerBlocks_;

  // Central FieldSet (for randomization)
  atlas::FieldSet centralFieldSet_;

  // Increment variables (to read ensemble members)
  oops::Variables incVars_;

  // Hybrid weights
  double scalarWeightSqrt_;
  atlas::FieldSet fileWeightSqrt_;
};

// -----------------------------------------------------------------------------

}  // namespace saber
