/*
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/detail/partitioner/TransPartitioner.h"

#include "oops/util/Logger.h"

#include "saber/interpolation/GaussToCS.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace saber {
namespace interpolation {

namespace {

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid) {
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(1));
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussToCS> makerGaussToCS_("gauss to cubed-sphere-dual");

// -----------------------------------------------------------------------------
// Note that this is slower than this needs to be
// as we are need to create 2 grid objects (very slow)
// In the future it might make sense to include an atlas grid (if available) from
// the model in outerGeometryData.
GaussToCS::GaussToCS(const oops::GeometryData & outerGeometryData,
                     const std::vector<size_t> & activeVariableSizes,
                     const oops::Variables & outerVars,
                     const Parameters_ & params,
                     const atlas::FieldSet & xb,
                     const atlas::FieldSet & fg,
                     const std::vector<atlas::FieldSet> & fsetVec)
  : innerVars_(outerVars),
    activeVars_(params.activeVariables.value().get_value_or(innerVars_)),
    CSFunctionSpace_(outerGeometryData.functionSpace()),
    gaussGrid_(params.gaussGridUid.value()),
    gaussFunctionSpace_(createGaussFunctionSpace(gaussGrid_)),
    gaussPartitioner_(new TransPartitioner()),
    csgrid_(CSFunctionSpace_.mesh().grid()),
    interp_(gaussPartitioner_, gaussFunctionSpace_, csgrid_, CSFunctionSpace_),
    innerGeometryData_(gaussFunctionSpace_, outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(),
                       outerGeometryData.comm())

{
  oops::Log::trace() << classname() << "::GaussToCS starting" << std::endl;
  oops::Log::trace() << classname() << "::GaussToCS done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
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

  // On input: fieldset on gaussian mesh

  // Create fieldset on cubed-sphere mesh.

  atlas::FieldSet csFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field csField =
      CSFunctionSpace_.createField<double>(
          atlas::option::name(fieldname) |
          atlas::option::levels(gaussFieldSet[fieldname].levels()) |
          atlas::option::halo(1));
    csFieldSet.add(csField);
  }

  // Interpolate to cubed sphere
  interp_.execute(gaussFieldSet, csFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    newFields.add(csFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiply done"
                     << fieldSet.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // On input: fieldset on gaussian grid
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet csFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       csFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // Create gauss fieldset
  atlas::FieldSet gaussFieldSet;
  for (const auto & fieldname : activeVars_.variables()) {
    atlas::Field gaussField =
      gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
          atlas::option::levels(csFieldSet[fieldname].levels()) |
          atlas::option::halo(1));
    gaussFieldSet.add(gaussField);
  }

  // Adjoint of interpolation from gauss to dual cubed sphere
  interp_.executeAdjoint(gaussFieldSet, csFieldSet);

  for (const auto & fieldname : activeVars_.variables()) {
    newFields.add(gaussFieldSet[fieldname]);
  }

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void GaussToCS::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
