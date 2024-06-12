/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/SuperMoistIncrOp.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SuperMoistIncrOp>
  makerSuperMoistIncrOp_("mo_super_mio");

// -----------------------------------------------------------------------------

SuperMoistIncrOp::SuperMoistIncrOp(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(getUnionOfInnerActiveAndOuterVars(params, outerVars)),
    intermediateTempVars_(params.intermediateTempVars(outerVars)),
    MIO_(std::make_unique<MoistIncrOp>(outerGeometryData,
                                       outerVars,
                                       covarConf,
                                       params.moistIncrOp,
                                       xb, fg)),
    exnerThetaToTemp_(std::make_unique<AirTemperature>(outerGeometryData,
                                                       MIO_->innerVars(),
                                                       covarConf,
                                                       params.airTemperature,
                                                       xb, fg))

{
  oops::Log::trace() << classname() << "::SuperMoistIncrOp starting" << std::endl;
  oops::Log::trace() << classname() << "::SuperMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

SuperMoistIncrOp::~SuperMoistIncrOp() {
  oops::Log::trace() << classname() << "::~SuperMoistIncrOp starting" << std::endl;
  util::Timer timer(classname(), "~SuperMoistIncrOp");
  oops::Log::trace() << classname() << "::~SuperMoistIncrOp done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  exnerThetaToTemp_->multiply(fset);
  MIO_->multiply(fset);
  fset.removeFields(intermediateTempVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  MIO_->multiplyAD(fset);
  exnerThetaToTemp_->multiplyAD(fset);
  fset.removeFields(intermediateTempVars_);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  exnerThetaToTemp_->multiply(fset);
  MIO_->leftInverseMultiply(fset);
  fset.removeFields(intermediateTempVars_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperMoistIncrOp::directCalibration(const oops::FieldSets & fsets) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  // Derive air temperature and total water perturbations first
  oops::FieldSets fsets_copy(fsets);
  for (size_t fset_index = 0; fset_index < fsets.ens_size(); fset_index++) {
    auto & fset = fsets_copy[fset_index];
    exnerThetaToTemp_->multiply(fset);
    MIO_->leftInverseMultiply(fset);
  }

  // Direct calibration of MIO block
  MIO_->directCalibration(fsets_copy);

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

void SuperMoistIncrOp::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
