/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/interpolation/Interpolation.h"

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
  std::vector<double> lats;
  std::vector<double> lons;
  geom.latlon(lats, lons, false);
  innerGeomData_->setGlobalTree(lats, lons);

  interp_.reset(new oops::GlobalInterpolator(
          params.localInterpConf.value(), *innerGeomData_,
          outerGeometryData.functionSpace(), outerGeometryData.comm()));

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.shape(1);
    atlas::Field field = outerGeomData_.functionSpace()->createField<double>(
        atlas::option::name(varName) | atlas::option::levels(nlev));
    field.metadata() = f.metadata();
    targetFieldSet.add(field);
  }

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    sourceFieldSet.add(fieldSet[varName]);
  }

  // Interpolate to target/outer grid
  interp_->apply(sourceFieldSet, targetFieldSet);

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

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.shape(1);
    atlas::Field field = innerGeomData_->functionSpace()->createField<double>(
        atlas::option::name(varName) | atlas::option::levels(nlev));
    field.metadata() = f.metadata();
    sourceFieldSet.add(field);
  }

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    targetFieldSet.add(fieldSet[varName]);
  }

  // Zero field
  util::zeroFieldSet(sourceFieldSet);

  // (Adjoint of:) Interpolate to target/outer grid
  interp_->applyAD(sourceFieldSet, targetFieldSet);

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
  if (!inverseInterp_) {
    inverseInterp_.reset(new oops::GlobalInterpolator(
          params_.localInterpConf.value(), outerGeomData_,
          innerGeomData_->functionSpace(), innerGeomData_->comm()));
  }

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.shape(1);
    atlas::Field field = innerGeomData_->functionSpace()->createField<double>(
        atlas::option::name(varName) | atlas::option::levels(nlev));
    field.metadata() = f.metadata();
    targetFieldSet.add(field);
  }

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    sourceFieldSet.add(fieldSet[varName]);
  }

  // Interpolate to target/inner grid
  inverseInterp_->apply(sourceFieldSet, targetFieldSet);

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
