/*
 * (C) Copyright 2022- UCAR
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralToGauss.h"

#include "oops/util/Logger.h"

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SpectralToGauss> makerSpectralToGauss_("spectral to gauss");

// -----------------------------------------------------------------------------

SpectralToGauss::SpectralToGauss(const oops::GeometryData & outerGeometryData,
                                 const std::vector<size_t> & activeVariableSizes,
                                 const oops::Variables & outerVars,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const std::vector<atlas::FieldSet> & fsetVec)
  : innerVars_(outerVars),
    activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    innerGeometryData_(specFunctionSpace_, outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(), outerGeometryData.comm())

{
  oops::Log::trace() << classname() << "::SpectralToGauss starting" << std::endl;
  oops::Log::trace() << classname() << "::SpectralToGauss done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet specFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       specFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // On input: fieldset in spectral space

  // Create fieldset on gaussian grid
  atlas::FieldSet gaussFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field gaussField =
      gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                 atlas::option::levels(specFieldSet[fieldname].levels()));
    gaussField.haloExchange();
    atlas::array::make_view<double, 2>(gaussField).assign(0.0);
    gaussFieldSet.add(gaussField);
  }

  // Transform to gaussian grid
  trans_.invtrans(specFieldSet, gaussFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    gaussFieldSet[fieldname].haloExchange();
    newFields.add(gaussFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << fieldSet.field_names() << std::endl;

  // On input: fieldset on gaussian grid
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       gaussFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // Create spectral fieldset
  atlas::FieldSet specFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field specField =
      specFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                 atlas::option::levels(gaussFieldSet[fieldname].levels()));
    specFieldSet.add(specField);
  }

  // Transform to spectral space
  trans_.invtrans_adj(gaussFieldSet, specFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    newFields.add(specFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
