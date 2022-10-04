/*
 * (C) Crown Copyright 2022 Met Office
 * (C) Copyright 2022- UCAR
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

SpectralToGauss::SpectralToGauss(const oops::GeometryData & outputGeometryData,
                                 const std::vector<size_t> & activeVariableSizes,
                                 const oops::Variables & outputVars,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const std::vector<atlas::FieldSet> & fsetVec)
  : gaussFunctionSpace_(outputGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    inputGeometryData_(specFunctionSpace_, outputGeometryData.fieldSet(),
                       outputGeometryData.levelsAreTopDown(), outputGeometryData.comm()),
    inputVars_(outputVars)
{
  oops::Log::trace() << classname() << "::SpectralToGauss starting" << std::endl;
  oops::Log::trace() << classname() << "::SpectralToGauss done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // On input: fieldset in spectral space
  const atlas::FieldSet & specFieldSet = fieldSet;

  // Create fieldset on gaussian grid
  atlas::FieldSet gaussFieldSet;
  for (const auto & fieldname : specFieldSet.field_names()) {
    atlas::Field gaussField =
      gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                 atlas::option::levels(specFieldSet[fieldname].levels()));
    gaussFieldSet.add(gaussField);
  }

  // Transform to gaussian grid
  trans_.invtrans(specFieldSet, gaussFieldSet);
  fieldSet = gaussFieldSet;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // On input: fieldset on gaussian grid
  const auto & gaussFieldSet = fieldSet;

  // Create spectral fieldset
  atlas::FieldSet specFieldSet;
  for (const auto & fieldname : gaussFieldSet.field_names()) {
    atlas::Field specField =
      specFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                 atlas::option::levels(gaussFieldSet[fieldname].levels()));
    specFieldSet.add(specField);
  }

  // Transform to spectral space
  trans_.invtrans_adj(gaussFieldSet, specFieldSet);
  fieldSet = specFieldSet;

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
