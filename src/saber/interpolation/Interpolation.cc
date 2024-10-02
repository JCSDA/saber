/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/Interpolation.h"

#include "atlas/util/Config.h"

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<Interpolation> makerInterpolation_("interpolation");

// -----------------------------------------------------------------------------

Interpolation::Interpolation(const oops::GeometryData & outerGeometryData,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const oops::FieldSet3D & xb,
                             const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    params_(params), outerGeomData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars))
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;

  // Set up GeometryData
  Geometry geom(params.innerGeom, outerGeometryData.comm());
  innerGeomData_.reset(new oops::GeometryData(geom.functionSpace(), geom.fields(),
                                              true, outerGeometryData.comm()));

  if (params.interpType.value() == "global") {
    globalInterp_.reset(new oops::GlobalInterpolator(
      params.forwardInterpConf.value(), *innerGeomData_,
      outerGeometryData.functionSpace(), outerGeometryData.comm()));
  } else if (params.interpType.value() == "regional") {
    regionalInterp_.reset(new atlas::Interpolation(
       atlas::util::Config("type", "regional-linear-2d"),
       innerGeomData_->functionSpace(), outerGeomData_.functionSpace()));
  } else {
    throw eckit::UserError("wrong interpolator type: " + params.interpType.value(), Here());
  }

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & var : activeVars_) {
    sourceFieldSet.add(fieldSet[var.name()]);
  }

  // Interpolate to target/outer grid
  atlas::FieldSet targetFieldSet;
  if (globalInterp_) {
    globalInterp_->apply(sourceFieldSet, targetFieldSet);
  }
  if (regionalInterp_) {
    for (const auto & var : activeVars_) {
      const atlas::Field sourceField = sourceFieldSet[var.name()];
      atlas::Field targetField = outerGeomData_.functionSpace().createField<double>(
          atlas::option::name(var.name()) | atlas::option::levels(sourceField.levels()));
      targetField.metadata() = sourceField.metadata();
      auto targetView = atlas::array::make_view<double, 2>(targetField);
      targetView.assign(0.0);
      targetFieldSet.add(targetField);
    }
    regionalInterp_->execute(sourceFieldSet, targetFieldSet);
  }

  // Add passive variables
  for (const auto & f : fieldSet) {
    if (!activeVars_.has(f.name())) {
      targetFieldSet.add(f);
    }
  }

  // Reset
  fieldSet.fieldSet() = targetFieldSet;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & var : activeVars_) {
    targetFieldSet.add(fieldSet[var.name()]);
  }

  // (Adjoint of:) Interpolate to target/outer grid
  atlas::FieldSet sourceFieldSet;
  if (globalInterp_) {
    globalInterp_->applyAD(sourceFieldSet, targetFieldSet);
  }
  if (regionalInterp_) {
    for (const auto & var : activeVars_) {
      const atlas::Field targetField = targetFieldSet[var.name()];
      atlas::Field sourceField = innerGeomData_->functionSpace().createField<double>(
          atlas::option::name(var.name()) | atlas::option::levels(targetField.levels()));
      sourceField.metadata() = targetField.metadata();
      auto sourceView = atlas::array::make_view<double, 2>(sourceField);
      sourceView.assign(0.0);
      sourceFieldSet.add(sourceField);
    }
    regionalInterp_->execute_adjoint(sourceFieldSet, targetFieldSet);
  }

  // Copy passive variables
  for (const auto & f : fieldSet) {
    if (!activeVars_.has(f.name())) {
      sourceFieldSet.add(f);
    }
  }

  fieldSet.fieldSet() = sourceFieldSet;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  // Prepare inverse interpolator
  if (!inverseGlobalInterp_ && globalInterp_) {
    inverseGlobalInterp_.reset(new oops::GlobalInterpolator(
      params_.inverseInterpConf.value(), outerGeomData_,
      innerGeomData_->functionSpace(), innerGeomData_->comm()));
  }
  if (!inverseRegionalInterp_ && regionalInterp_) {
    inverseRegionalInterp_.reset(new atlas::Interpolation(
       atlas::util::Config("type", "regional-linear-2d"),
       outerGeomData_.functionSpace(), innerGeomData_->functionSpace()));
  }

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & var : activeVars_) {
    sourceFieldSet.add(fieldSet[var.name()]);
  }

  // Interpolate to target/inner grid
  atlas::FieldSet targetFieldSet;
  if (inverseGlobalInterp_) {
    inverseGlobalInterp_->apply(sourceFieldSet, targetFieldSet);
  }
  if (inverseRegionalInterp_) {
    for (const auto & var : activeVars_) {
      const atlas::Field sourceField = sourceFieldSet[var.name()];
      atlas::Field targetField = outerGeomData_.functionSpace().createField<double>(
          atlas::option::name(var.name()) | atlas::option::levels(sourceField.levels()));
      targetField.metadata() = sourceField.metadata();
      auto targetView = atlas::array::make_view<double, 2>(targetField);
      targetView.assign(0.0);
      targetFieldSet.add(targetField);
    }
    inverseRegionalInterp_->execute(sourceFieldSet, targetFieldSet);
  }

  // Add passive variables
  for (const auto & f : fieldSet) {
    if (!activeVars_.has(f.name())) {
      targetFieldSet.add(f);
    }
  }

  // Reset
  fieldSet.fieldSet() = targetFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D Interpolation::generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                                      const oops::Variables & innerVars) const {
  oops::FieldSet3D fset(this->validTime(), innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D Interpolation::generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                                      const oops::Variables & outerVars) const {
  oops::FieldSet3D fset(this->validTime(), outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
