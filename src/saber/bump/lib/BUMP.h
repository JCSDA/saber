/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/log/Channel.h"

#include "saber/bump/lib/type_bump.h"

namespace bump_lib {

// -----------------------------------------------------------------------------

class BUMP {
 public:
  // Constructor
  BUMP(const eckit::mpi::Comm &,
       eckit::Channel &,
       eckit::Channel &,
       const atlas::FunctionSpace &,
       const atlas::FieldSet &,
       const std::vector<size_t> &,
       const std::vector<std::string> &,
       const eckit::Configuration &,
       const eckit::Configuration &);

  // Destructor
  ~BUMP();

  // Read input ATLAS files
  void readAtlasFiles();

  // Fields to read
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead(
    const std::vector<eckit::LocalConfiguration>);

  // Return input pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> inputs() const
    {return inputs_;}

  // Add fields
  void addField(const atlas::FieldSet &);

  // Add ensemble
  void addEnsemble(const std::vector<atlas::FieldSet> &);

  // Dual resolution setup
  void dualResolutionSetup(const atlas::FunctionSpace &,
                           const atlas::FieldSet &);

  // Iterative update
  void iterativeUpdate(const atlas::FieldSet &, const size_t &);

  // Write output ATLAS files
  void writeAtlasFiles() const;

  // Fields to write
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite(
    const std::vector<eckit::LocalConfiguration>) const;

  // Fortran interfaces
  void runDrivers();
  void multiplyVbal(atlas::FieldSet &) const;
  void multiplyVbalAd(atlas::FieldSet &) const;
  void inverseMultiplyVbal(atlas::FieldSet &) const;
  void multiplyStdDev(atlas::FieldSet &) const;
  void inverseMultiplyStdDev(atlas::FieldSet &) const;
  void randomizeNicas(atlas::FieldSet &) const;
  void multiplyNicas(atlas::FieldSet &) const;
  void multiplyPsiChiToUV(atlas::FieldSet &) const;
  void multiplyPsiChiToUVAd(atlas::FieldSet &) const;
  size_t getCvSize() const;
  void multiplyNicasSqrt(const atlas::Field &, atlas::FieldSet &, const size_t &) const;
  void multiplyNicasSqrtAd(const atlas::FieldSet &, atlas::Field &, const size_t &) const;

 private:
  std::vector<int> keyBUMP_;
  const eckit::mpi::Comm * comm_;
  eckit::Channel * infoChannel_;
  eckit::Channel * testChannel_;
  atlas::FunctionSpace fspace_;
  std::vector<size_t> variableSizes_;
  std::vector<std::string> vars_;
  eckit::LocalConfiguration covarConf_;
  eckit::LocalConfiguration bumpConf_;
  std::vector<size_t> nens_;
  bool iterativeEnsembleLoading_;
  bool waitForDualResolution_;
  std::string gridUid_;
  std::string dualResolutionGridUid_;
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> inputs_;

  eckit::LocalConfiguration getFileConf(const eckit::mpi::Comm &,
                                        const eckit::Configuration &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bump_lib
