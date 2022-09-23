/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SPTOUVParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SPTOUVParameters, SaberOuterBlockParametersBase)
 public:
  oops::RequiredParameter<std::string> gaussGridUid{"gauss_grid_uid",
    "Gauss Grid UID", this};
};

// -----------------------------------------------------------------------------


class SPTOUV : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::spectralb::SPTOUV";}

  typedef SPTOUVParameters Parameters_;

  SPTOUV(const eckit::mpi::Comm & comm,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<std::size_t> &,
         const eckit::Configuration & conf,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);

  virtual ~SPTOUV();

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  atlas::FunctionSpace outputFunctionSpace_;
  atlas::StructuredGrid gaussGrid_;
  oops::Variables inputVars_;
  oops::Variables activeVars_;
  std::vector<std::size_t> activeVariableSizes_;
  atlas::functionspace::Spectral specFS_;
  atlas::trans::Trans transFS_;

  void print(std::ostream &) const override;

};

}  // namespace spectralb
}  // namespace saber
