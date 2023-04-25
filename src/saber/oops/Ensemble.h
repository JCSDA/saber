/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"

#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class EnsembleParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(EnsembleParameters, SaberBlockParametersBase)

 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class Ensemble : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::Ensemble";}

  typedef EnsembleParameters Parameters_;

  Ensemble(const oops::GeometryData &,
           const std::vector<size_t> &,
           const oops::Variables &,
           const eckit::Configuration &,
           const Parameters_ &,
           const atlas::FieldSet &,
           const atlas::FieldSet &,
           const size_t &);

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;

  void setLocalization(std::unique_ptr<SaberBlockChain>) override;

 private:
  std::vector<atlas::FieldSet> ensemble_;
  std::unique_ptr<SaberBlockChain> loc_;
  size_t timeRank_;
  const oops::Variables vars_;
  const eckit::mpi::Comm & comm_;
  int seed_ = 7;  // For reproducibility
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
