/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <map>

#include "saber/interpolation/VertProj.h"

#include "eckit/config/LocalConfiguration.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace interpolation {

namespace {

oops::Variables createInnerVars(
    const atlas::idx_t & innerVerticalLevels,
    const oops::Variables & activeVars,
    const oops::Variables & outerVars) {
  oops::Variables innerVars(outerVars);
  for (const auto & var : activeVars) {
    if (!innerVars.has(var)) {
      innerVars.push_back(var);
    }
    innerVars[var.name()].setLevels(innerVerticalLevels);
  }
  return innerVars;
}

void verticalProjection(const atlas::FieldSet & fsetin, atlas::FieldSet & fsetout) {
  for (const auto & f : fsetin) {
    auto fldInView = atlas::array::make_view<const double, 2>(f);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetout[f.name()]);
    for (atlas::idx_t jn = 0; jn < fldOutView.shape()[0]; ++jn) {
      for (atlas::idx_t jl = 0; jl < fldOutView.shape()[1]; ++jl) {
        const atlas::idx_t jl2 = static_cast<atlas::idx_t>(jl * fldInView.shape()[1]
                                 / fldOutView.shape()[1]);
        fldOutView(jn, jl) = fldInView(jn, jl2);
      }
    }
  }
}

void verticalProjectionInverse(atlas::FieldSet & fsetin, const atlas::FieldSet & fsetout) {
  for (auto & f : fsetin) {
    auto fldInView = atlas::array::make_view<double, 2>(f);
    auto fldOutView = atlas::array::make_view<const double, 2>(fsetout[f.name()]);
    for (atlas::idx_t jn = 0; jn < fldOutView.shape()[0]; ++jn) {
      for (atlas::idx_t jl = 0; jl < fldOutView.shape()[1]; ++jl) {
        const atlas::idx_t jl2 = static_cast<atlas::idx_t>(jl * fldInView.shape()[1]
                                 / fldOutView.shape()[1]);
        fldInView(jn, jl2) = fldOutView(jn, jl);
      }
    }
  }
}

void verticalProjectionAD(atlas::FieldSet & fsetin, atlas::FieldSet & fsetout) {
  for (auto & f : fsetin) {
    auto fldInView = atlas::array::make_view<double, 2>(f);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetout[f.name()]);
    for (atlas::idx_t jn = 0; jn < fldOutView.shape()[0]; ++jn) {
      for (atlas::idx_t jl = fldOutView.shape()[1]-1; jl >= 0; --jl) {
        atlas::idx_t jl2 = jl * fldInView.shape()[1] /fldOutView.shape()[1];
        fldInView(jn, jl2) += fldOutView(jn, jl);
        fldOutView(jn, jl) = 0.0;
      }
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VertProj> makerVertProj_("simple vertical projection");

// -----------------------------------------------------------------------------
// Note that this is assumes that the variables in activevars
//      exist in the variable level mappings of the inner and outer
//      GeometryData objects.
VertProj::VertProj(const oops::GeometryData & outerGeometryData,
                   const oops::Variables & outerVars,
                   const eckit::Configuration & covarConf,
                   const Parameters_ & params,
                   const oops::FieldSet3D & xb,
                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    outerGeometryData_(outerGeometryData),
    outerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    innerVars_(createInnerVars(params.innerVerticalLevels, activeVars_, outerVars))
{
  oops::Log::trace() << classname() << "::VertProj starting" << std::endl;
  oops::Log::trace() << classname() << "::VertProj done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertProj::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet vertFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       vertFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // Create fieldset on functionspace
  atlas::FieldSet modelFieldSet;
  // Note: ambiguity in levels in active variables; get levels from outer vars
  for (const auto & var : activeVars_) {
    atlas::Field modelField =
      outerGeometryData_.functionSpace().createField<double>(
          atlas::option::name(var.name()) |
          atlas::option::levels(outerVars_[var.name()].getLevels()) |
          atlas::option::halo(1));
    atlas::array::make_view<double, 2>(modelField).assign(0.0);
    modelField.set_dirty(false);
    modelFieldSet.add(modelField);
  }

  // Simple prologation scheme
  verticalProjection(vertFieldSet, modelFieldSet);

  for (const auto & var : activeVars_) {
    newFields.add(modelFieldSet[var.name()]);
  }

  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname() << "::multiply done"
                     << fieldSet.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void VertProj::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // On input: fieldset on model levels
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet modelFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       modelFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // Create vert fieldset
  atlas::FieldSet vertFieldSet;
  // Note: ambiguity in levels in active variables; get levels from inner vars
  for (const auto & var : activeVars_) {
    atlas::Field vertField =
      outerGeometryData_.functionSpace().createField<double>
        (atlas::option::name(var.name()) |
         atlas::option::levels(innerVars_[var.name()].getLevels()) |
         atlas::option::halo(1));
    atlas::array::make_view<double, 2>(vertField).assign(0.0);
    vertField.set_dirty(false);
    vertFieldSet.add(vertField);
  }

  // Adjoint of simple prolongation scheme
  verticalProjectionAD(vertFieldSet, modelFieldSet);

  for (const auto & var : activeVars_) {
    newFields.add(vertFieldSet[var.name()]);
  }

  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname()
                     << "::multiplyAD done" << fieldSet.field_names()  << std::endl;
}


// -----------------------------------------------------------------------------

void VertProj::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname()
                     << "::leftInverseStarting starting" << std::endl;

  // On input: fieldset on model levels
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet modelFieldSet = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fieldSet.field_names();

  for (auto & fieldname : fieldSet.field_names()) {
     if (activeVars_.has(fieldname)) {
       modelFieldSet.add(fieldSet[fieldname]);
     } else {
       newFields.add(fieldSet[fieldname]);
     }
  }

  // Create vert fieldset
  atlas::FieldSet vertFieldSet;
  for (const auto & var : activeVars_) {
    atlas::Field vertField =
      outerGeometryData_.functionSpace().createField<double>
        (atlas::option::name(var.name()) |
         atlas::option::levels(innerVars_[var.name()].getLevels()) |
         atlas::option::halo(1));
    atlas::array::make_view<double, 2>(vertField).assign(0.0);
    vertField.set_dirty(false);
    vertFieldSet.add(vertField);
  }

  // left inverse mutiply
  verticalProjectionInverse(vertFieldSet, modelFieldSet);

  for (const auto & var : activeVars_) {
    newFields.add(vertFieldSet[var.name()]);
  }

  fieldSet.fieldSet() = newFields;

  oops::Log::trace() << classname()
                     << "::leftInverseMultiply done"
                     << fieldSet.field_names()  << std::endl;
}

// -----------------------------------------------------------------------------

void VertProj::print(std::ostream & os) const {
  os << classname();
}

}  // namespace interpolation
}  // namespace saber
