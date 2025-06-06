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
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

namespace {

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
    auto innerView = atlas::array::make_view<double, 2>(fsetInput[key]);
    for (const auto & component : v) {
      auto outerView = atlas::array::make_view<double, 2>(fsetOutput[component.name()]);
      if (levelsAreTopDown && outerView.shape(1) == 1) {
        // jn is horizontal index
        for (atlas::idx_t jn = 0; jn < outerView.shape(0); ++jn) {
          // copy FROM last index of inner field
          outerView(jn, 0) = innerView(jn, innerView.shape(1)-1);
        }
      } else {
        // jn is horizontal index
        for (atlas::idx_t jn = 0; jn < outerView.shape(0); ++jn) {
          // loop over vertical levels that exist in OUTER field
          for (atlas::idx_t jl = 0; jl < outerView.shape(1); ++jl) {
            // copy in place
            outerView(jn, jl) = innerView(jn, jl);
          }
        }
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
    auto innerView = atlas::array::make_view<double, 2>(fsetOutput[key]);
    for (const auto & component : v) {
      if (!interpType.empty()) {
        // check other fields have same 'interp_type' as component0
        ASSERT(fsetInput[component.name()].metadata().has("interp_type"));
        ASSERT(fsetInput[component.name()].metadata().get<std::string>("interp_type")
                                                                      == interpType);
      }
      auto outerView = atlas::array::make_view<double, 2>(fsetInput[component.name()]);
      // assert same horizontal size for inner and outer fields
      ASSERT(outerView.shape(0) == innerView.shape(0));
      // assert same vertical size OR outer Field has one level (i.e., is a surface field)
      ASSERT(outerView.shape(1) == innerView.shape(1) || outerView.shape(1) == 1);
      if (levelsAreTopDown && outerView.shape(1) == 1) {
        // candidate for openMP parallel for??
        for (atlas::idx_t jn = 0; jn < outerView.shape(0); ++jn) {
          // gather to end of array
          innerView(jn, innerView.shape(1)-1) += outerView(jn, 0);
          outerView(jn, 0) = 0.0;
        }
      } else {
        // candidate for openMP parallel for??
        for (atlas::idx_t jn = 0; jn < outerView.shape(0); ++jn) {
          // only loop over levels that exist in incoming field
          for (atlas::idx_t jl = 0; jl < outerView.shape(1); ++jl) {
            innerView(jn, jl) += outerView(jn, jl);
            outerView(jn, jl) = 0.0;
          }
        }
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
