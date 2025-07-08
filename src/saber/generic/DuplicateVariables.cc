/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/DuplicateVariables.h"

#include <algorithm>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

namespace {

using View = atlas::array::LocalView<double, 1>;
using ConstView = atlas::array::LocalView<const double, 1>;

oops::Variables createActiveVars(const std::vector<VariableGroupParameters> & gps,
                                 const oops::Variables & outerVars)  {
  oops::Variables activeVars;
  for (const VariableGroupParameters & gp : gps) {
    // v is read from yaml, and does not have levels metadata; copy variables
    // from outerVars that have relevant metadata.
    oops::Variables v = gp.groupComponents.value();
    ASSERT(v.size() > 0);
    std::vector<int> lvls;
    for (const auto & var : v) {
      activeVars.push_back(outerVars[var.name()]);
      lvls.push_back(outerVars[var.name()].getLevels());
    }
    auto max_level = std::max_element(lvls.begin(), lvls.end());
    // add the group variable name and assign the same levels as in other vars
    std::string key = gp.groupVariableName.value();
    eckit::LocalConfiguration conf;
    conf.set("levels", *max_level);
    activeVars.push_back({key, conf});
  }

  return activeVars;
}


oops::Variables createInnerVars(const oops::Variables & outerVars,
                                const oops::Variables & activeVars,
                                const std::vector<VariableGroupParameters> & gps) {
  oops::Variables innerVars;
  for (const auto & var : outerVars) {
    if (!activeVars.has(var)) {
      innerVars.push_back(var);
    }
  }

  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    innerVars.push_back(activeVars[key]);
  }

  return innerVars;
}
/**
* @brief Broadcasts (copies) the grouped field from the
  *      INNER fieldset to all individual fields in group.
  *      Used in forward multiply() which takes INNER -> OUTER
  *      (assumes all activeFields are allocated)
  *
  * @param gps vector of configurations prescribing fields to be grouped
  * @param fsetInput the input fieldset of saber INNER variables
  * @param fsetOutput the output fields set with saber OUTER variable
  */
void copyFields(const std::vector<VariableGroupParameters> & gps,
                const atlas::FieldSet & fsetInput,
                atlas::FieldSet & fsetOutput,
                bool levelsAreTopDown) {
  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    oops::Variables v = gp.groupComponents.value();
    // copy metadata in to re-expanded fields
    std::string interpType;
    if (fsetInput[key].metadata().has("interp_type")) {
      interpType = fsetInput[key].metadata().get<std::string>("interp_type");
    }
    for (const auto & component : v) {
      if (levelsAreTopDown && fsetOutput[component.name()].shape(1) == 1) {
        // surface variable is at end of array
        util::for_each_column(
          util::IndexRange::include_halo,  // TODO(NC): when upgrade to atlas0.43 delete this line
          // copy surface FROM last index of inner field
          [](ConstView inner, View outer) {
            const atlas::idx_t surfLevel = inner.shape(0)-1;
            outer(0) = inner(surfLevel);
          },
          fsetInput[key], fsetOutput[component.name()]);
      } else {
        // only copy values up to max level in OUTER field
        util::for_each_column(
          util::IndexRange::include_halo,  // TODO(NC): when upgrade to atlas0.43 delete this line
          [](ConstView inner, View outer) {
            for (atlas::idx_t jl = 0; jl < outer.shape(0); ++jl){ outer(jl) = inner(jl); }
          },
            fsetInput[key], fsetOutput[component.name()]);
      }
      if (!interpType.empty()) {
        fsetOutput[component.name()].metadata().set("interp_type", interpType);
      }
    }
  }
}

/**
 * @brief Gathers (+='s) prescribed groups of fields into a
 *        single field with name set by the string 'key'.
 *        Used in adjoint direction, i.e. multiplyAD()
 *        which takes OUTER -> INNER
 *
 * @param gps vector of configurations prescribing fields to be grouped
 * @param fsetInput the input fieldset of saber OUTER variables
 * @param fsetOutput the output fields set with saber INNER variable
 */
void gatherFields(const std::vector<VariableGroupParameters> & gps,
                  atlas::FieldSet & fsetInput,
                  atlas::FieldSet & fsetOutput,
                  bool levelsAreTopDown) {
  for (const VariableGroupParameters & gp : gps) {
    std::string key = gp.groupVariableName.value();
    oops::Variables v = gp.groupComponents.value();
    std::string interpType;

    // copy metadata from 0th variable, if present
    if (fsetInput[v[0].name()].metadata().has("interp_type")) {
      interpType = fsetInput[v[0].name()].metadata().get<std::string>("interp_type");
    }
    for (const auto & component : v) {
      if (!interpType.empty()) {
        // check other fields have same 'interp_type' as component0
        ASSERT(fsetInput[component.name()].metadata().has("interp_type"));
        ASSERT(fsetInput[component.name()].metadata().get<std::string>("interp_type")
                                                                      == interpType);
      }
      // assert same horizontal size for inner and outer fields
      ASSERT(fsetInput[component.name()].shape(0) == fsetOutput[key].shape(0));
      // assert same vertical size OR outer Field has one level (i.e., is a surface field)
      ASSERT(fsetInput[component.name()].shape(1) == fsetOutput[key].shape(1)
             || fsetInput[component.name()].shape(1) == 1);
      if (levelsAreTopDown && fsetInput[component.name()].shape(1) == 1) {
        // surface variable is at end of array
        util::for_each_column(
          util::IndexRange::include_halo,  // TODO(NC): when upgrade to atlas0.43 delete this line
          // gather to end of array
          [](View inner, View outer) {
            inner(inner.shape(0)-1) += outer(0);
            outer.assign(0);
          },
          fsetOutput[key], fsetInput[component.name()]);
      } else {
        using View = atlas::array::LocalView<double, 1>;
        util::for_each_column(
          util::IndexRange::include_halo,  // TODO(NC): when upgrade to atlas0.43 delete this line
          [](View inner, View outer) {
            // only gather values from OUTER field up to max level in Input/Outer field
            for (atlas::idx_t jl = 0; jl < outer.shape(0); ++jl){
              inner(jl) += outer(jl);
              outer(jl) = 0; }
          },
          fsetOutput[key], fsetInput[component.name()]);
      }
    }
    if (!interpType.empty()) {
      fsetOutput[key].metadata().set("interp_type", interpType);
    }
  }
}

}  //  namespace


// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<DuplicateVariables>
    makerDuplicateVariables_("duplicate variables");

// -----------------------------------------------------------------------------

DuplicateVariables::DuplicateVariables(const oops::GeometryData & outerGeometryData,
                                       const oops::Variables & outerVars,
                                       const eckit::Configuration & covarConfig,
                                       const Parameters_ & params,
                                       const oops::FieldSet3D & xb,
                                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    groups_(params.variableGroupParameters.value()),
    outerVars_(outerVars),
    activeVars_(createActiveVars(groups_, outerVars_)),
    innerVars_(createInnerVars(outerVars_, activeVars_, groups_)),
    innerGeometryData_(outerGeometryData)
{
  oops::Log::trace() << classname() << "::DuplicateVariables starting" << std::endl;
  oops::Log::trace() << classname() << "::DuplicateVariables done" << std::endl;
}

// -----------------------------------------------------------------------------
void DuplicateVariables::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // The aim is to copy groupVars_ fields into component fields.
  // Some manipulation of fieldsets is needed for this.
  atlas::FieldSet fsetOut = atlas::FieldSet();

  for (const VariableGroupParameters & gp : groups_) {
    oops::Variables v = gp.groupComponents.value();
    for (const auto & component : v) {
      const size_t nlev = activeVars_[component.name()].getLevels();
      atlas::Field field =
        innerGeometryData_.functionSpace()->createField<double>(
        atlas::option::name(component.name()) | atlas::option::levels(nlev));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      field.set_dirty(false);
      fsetOut.add(field);
    }
  }

  // copy group fields into component fields
  copyFields(groups_, fset.fieldSet(), fsetOut, innerGeometryData_.levelsAreTopDown());

  // keep passive vars
  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateVariables::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // The aim is to do the adjoint of the copy groupVars_ fields into component fields.
  // This involves summing at the component fields together in
  // Some manipulation of fieldsets is needed for this.
  atlas::FieldSet fsetOut = atlas::FieldSet();

  // allocate group vars.
  for (const VariableGroupParameters & gp : groups_) {
    std::string key = gp.groupVariableName.value();
    const size_t nlev = activeVars_[key].getLevels();
    atlas::Field field =
      innerGeometryData_.functionSpace()->createField<double>(
        atlas::option::name(key) | atlas::option::levels(nlev));
    atlas::array::make_view<double, 2>(field).assign(0.0);
    field.set_dirty(false);
    fsetOut.add(field);
  }

  // sum component fields into group fields
  gatherFields(groups_, fset.fieldSet(), fsetOut, innerGeometryData_.levelsAreTopDown());

  // keep passive vars
  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiplyAD done = " << std::endl;
}

// -----------------------------------------------------------------------------

void DuplicateVariables::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
