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
                             const std::vector<size_t> & activeVariableSizes,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    params_(params), outerGeomData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars))
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;

  // Set up GeometryData
  //
  // We do this without going through the oops::Geometry<MODEL>, for two reasons:
  // 1. It avoids needing to make a dummy MODEL trait just to use Atlas. (Note: we then essentially
  //    negate this point by using quench "model" code to read a yaml into Atlas objects...)
  // 2. It works around a silly incompatibility between Atlas's gaussian grids and the KDTree data-
  //    structures used in the GeometryData. Specifically, Atlas gaussian grids place halo points
  //    "across the pole", i.e., with latitude coordinates ranging outside [-90,90], and this breaks
  //    the coordinate transforms use by the KDTree. Note this incompatibility affects only the
  //    local tree in the GeometryData, because that's the tree constructed using halo points. And
  //    since we don't actually need the local tree for the AtlasInterpolator used in this saber
  //    block, we can simply skip the initializiation of this data structure and conveniently avoid
  //    the incompatibility... an absolute house of cards!
  Geometry geom(params.innerGeom, outerGeometryData.comm());
  innerGeomData_.reset(new oops::GeometryData(geom.functionSpace(), outerGeometryData.fieldSet(),
                                              true, outerGeometryData.comm()));
  std::vector<double> lats;
  std::vector<double> lons;
  geom.latlon(lats, lons, false);
  innerGeomData_->setGlobalTree(lats, lons);

  interp_.reset(new oops::GlobalAtlasInterpolator(
        params.localInterpConf.value(), *innerGeomData_,
        outerGeometryData.functionSpace(), outerGeometryData.comm()));

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.levels();
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
  fieldSet = targetFieldSet;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Temporary FieldSet of active variables for interpolation source
  atlas::FieldSet sourceFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.levels();
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

  fieldSet = sourceFieldSet;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::leftInverseMultiply(atlas::FieldSet & fieldSet) const {
  if (!inverseInterp_) {
    inverseInterp_.reset(new oops::GlobalAtlasInterpolator(
          params_.localInterpConf.value(), outerGeomData_,
          innerGeomData_->functionSpace(), innerGeomData_->comm()));
  }

  // Temporary FieldSet of active variables for interpolation target
  atlas::FieldSet targetFieldSet;
  for (const auto & varName : activeVars_.variables()) {
    const auto & f = fieldSet[varName];
    const size_t nlev = f.levels();
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
  fieldSet = targetFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
