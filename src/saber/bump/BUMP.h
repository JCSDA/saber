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

#include "oops/base/FieldSets.h"
#include "oops/base/GeometryData.h"

#include "saber/bump/BUMPParameters.h"
#include "saber/bump/type_bump.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

class BUMP {
 public:
  static const std::string classname() {return "saber::bump::BUMP";}

  // Constructor
  BUMP(const oops::GeometryData &,
       const oops::Variables &,
       const eckit::Configuration &,
       const BUMPParameters &,
       const eckit::LocalConfiguration &,
       const oops::FieldSet3D &);

  // Destructor
  ~BUMP();

  // Read input ATLAS files
  void readAtlasFiles();

  // Fields to read
  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs(
    const std::vector<eckit::LocalConfiguration>);

  // Add fields
  void addField(const oops::FieldSet3D &);

  // Add ensemble
  void addEnsemble(const oops::FieldSets &);

  // Dual resolution setup
  void dualResolutionSetup(const atlas::FunctionSpace &,
                           const atlas::FieldSet &);

  // Iterative update
  void iterativeUpdate(const oops::FieldSet3D &, const size_t &);

  // Write output ATLAS files
  void writeAtlasFiles() const;

  // Fields to write
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite(
    const std::vector<eckit::LocalConfiguration>) const;

  // Fortran interfaces
  void runDrivers();
  void multiplyVbal(oops::FieldSet3D &) const;
  void multiplyVbalAd(oops::FieldSet3D &) const;
  void inverseMultiplyVbal(oops::FieldSet3D &) const;
  void multiplyStdDev(oops::FieldSet3D &) const;
  void inverseMultiplyStdDev(oops::FieldSet3D &) const;
  void randomizeNicas(oops::FieldSet3D &) const;
  void multiplyNicas(oops::FieldSet3D &) const;
  void multiplyPsiChiToUV(oops::FieldSet3D &) const;
  void multiplyPsiChiToUVAd(oops::FieldSet3D &) const;
  size_t getCvSize() const;
  void multiplyNicasSqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const;
  void multiplyNicasSqrtAd(const oops::FieldSet3D &, atlas::Field &, const size_t &) const;

 private:
  std::vector<int> keyBUMP_;
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace fspace_;
  const oops::Variables vars_;
  const util::DateTime validTime_;
  eckit::LocalConfiguration covarConf_;
  eckit::LocalConfiguration bumpConf_;
  std::vector<size_t> nens_;
  bool iterativeEnsembleLoading_;
  bool waitForDualResolution_;
  std::string gridUid_;
  std::string dualResolutionGridUid_;

  eckit::LocalConfiguration getFileConf(const eckit::mpi::Comm &,
                                        const eckit::Configuration &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
