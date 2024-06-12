/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <string>

#include "saber/generic/VertLocInterp.h"

#include "eckit/config/LocalConfiguration.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace vader {

namespace {


oops::Variables createInnerVars(
    const atlas::idx_t & innerVerticalLevels,
    const oops::Variables & activeVars,
    const oops::Variables & outerVars) {
  oops::Variables innerVars(outerVars);
  for (auto & var : activeVars) {
    innerVars[var.name()].setLevels(innerVerticalLevels);
  }
  return innerVars;
}


void verticalInterpolationFromFullLevels(const oops::Variables & activeVars,
                                         const atlas::FieldSet & fsetIn,
                                         atlas::FieldSet & fsetOut) {
  // winds are on half-levels
  oops::Variables windnames(std::vector<std::string>{"eastward_wind", "northward_wind"});

  // inner stagger is full-levels
  for (const auto & v : activeVars) {
    auto fldInView = atlas::array::make_view<const double, 2>(fsetIn[v.name()]);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetOut[v.name()]);

    if (((v.name().find("_levels_minus_one")) != std::string::npos) ||
        windnames.has(v)) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        fldOutView(jn, 0) = fldInView(jn, 0);
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]-1; ++jl) {
          fldOutView(jn, jl+1) = 0.5 * (fldInView(jn, jl) + fldInView(jn, jl+1));
        }
      }
    } else if ((v.name().find("_levels")) != std::string::npos) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        fldOutView(jn, 0) = fldInView(jn, 0);
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]-1; ++jl) {
          fldOutView(jn, jl+1) = 0.5 * (fldInView(jn, jl) + fldInView(jn, jl+1));
        }
        fldOutView(jn, fldOutView.shape()[1]-1) =
          fldInView(jn, fldInView.shape()[1]-1);
      }
    } else {
      // full levels
      fldOutView.assign(fldInView);
    }
  }
}

// TO DO(Marek) - write non-buggy version using height or eta coordinates.
void mo_buggy_verticalInterpolationFromHalfLevels(const oops::Variables & activeVars,
                                                  const atlas::FieldSet & fsetIn,
                                                  atlas::FieldSet & fsetOut) {
  // winds are on half-levels
  oops::Variables windnames(std::vector<std::string>{"eastward_wind", "northward_wind"});

  // The Met Office bug unfortunately uses the interpolation weights needed for
  // going from full levels to half levels, when in fact we are wanting to interpolate
  // from half-levels.  This is not an issue if the vertical staggering is uniform.
  // However in practice it is not.
  for (const auto & v : activeVars) {
    auto fldInView = atlas::array::make_view<const double, 2>(fsetIn[v.name()]);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetOut[v.name()]);

    if (((v.name().find("_levels_minus_one")) != std::string::npos) ||
        windnames.has(v)) {
      fldOutView.assign(fldInView);
    } else if ((v.name().find("_levels")) != std::string::npos) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]; ++jl) {
          fldOutView(jn, jl) = fldInView(jn, jl);
        }
        fldOutView(jn, fldOutView.shape()[1]-1) =
          fldInView(jn, fldInView.shape()[1]-1);
      }
    } else {
      // full levels
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]-1; ++jl) {
          fldOutView(jn, jl) = 0.5 * (fldInView(jn, jl) + fldInView(jn, jl+1));
        }
        fldOutView(jn, fldOutView.shape()[1]-1) =
          fldInView(jn, fldInView.shape()[1]-1);
      }
    }
  }
}


void verticalInterpolationFromFullLevelsAD(const oops::Variables & activeVars,
                                           atlas::FieldSet & fsetIn,
                                           atlas::FieldSet & fsetOut) {
  // winds are on half-levels
  oops::Variables windnames(std::vector<std::string>{"eastward_wind", "northward_wind"});

  for (const auto & v : activeVars) {
    auto fldInView = atlas::array::make_view<double, 2>(fsetIn[v.name()]);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetOut[v.name()]);

    if ((v.name().find("_levels_minus_one")) != std::string::npos ||
        windnames.has(v)) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        for (atlas::idx_t jl = fldInView.shape()[1]-2; jl > -1; --jl) {
          fldInView(jn, jl+1) += 0.5 * fldOutView(jn, jl+1);
          fldInView(jn, jl) += 0.5 * fldOutView(jn, jl+1);
          fldOutView(jn, jl+1) = 0.0;
        }
        fldInView(jn, 0) += fldOutView(jn, 0);
        fldOutView(jn, 0) = 0.0;
      }
    } else if ((v.name().find("_levels")) != std::string::npos) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        fldInView(jn, fldInView.shape()[1]-1) +=
          fldOutView(jn, fldOutView.shape()[1]-1);
        fldOutView(jn, fldOutView.shape()[1]-1) = 0.0;
        for (atlas::idx_t jl = fldInView.shape()[1]-2; jl > -1; --jl) {
          fldInView(jn, jl+1) += 0.5 * fldOutView(jn, jl+1);
          fldInView(jn, jl) += 0.5 * fldOutView(jn, jl+1);
          fldOutView(jn, jl+1) = 0.0;
        }
        fldInView(jn, 0) += fldOutView(jn, 0);
        fldOutView(jn, 0) = 0.0;
      }
    } else {
      // full levels
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]; ++jl) {
          fldInView(jn, jl) += fldOutView(jn, jl);
          fldOutView(jn, jl) = 0.0;
        }
      }
    }
  }
}


void mo_buggy_verticalInterpolationFromHalfLevelsAD(const oops::Variables & activeVars,
                                                    atlas::FieldSet & fsetIn,
                                                    atlas::FieldSet & fsetOut) {
  // winds are on half-levels
  oops::Variables windnames(std::vector<std::string>{"eastward_wind", "northward_wind"});

  for (const auto & v : activeVars) {
    auto fldInView = atlas::array::make_view<double, 2>(fsetIn[v.name()]);
    auto fldOutView = atlas::array::make_view<double, 2>(fsetOut[v.name()]);

    if ((v.name().find("_levels_minus_one")) != std::string::npos ||
        windnames.has(v)) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]; ++jl) {
          fldInView(jn, jl) += fldOutView(jn, jl);
          fldOutView(jn, jl) = 0.0;
        }
      }
    } else if ((v.name().find("_levels")) != std::string::npos) {
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        fldInView(jn, fldInView.shape()[1]-1) +=
          fldOutView(jn, fldOutView.shape()[1]-1);
        fldOutView(jn, fldOutView.shape()[1]-1) = 0.0;
        for (atlas::idx_t jl = 0; jl < fldInView.shape()[1]; ++jl) {
          fldInView(jn, jl) += fldOutView(jn, jl);
          fldOutView(jn, jl) = 0.0;
        }
      }
    } else {
      // full levels
      for (atlas::idx_t jn = 0; jn < fldInView.shape()[0]; ++jn) {
        fldInView(jn, fldInView.shape()[1]-1) +=
          fldOutView(jn, fldOutView.shape()[1]-1);
        fldOutView(jn, fldOutView.shape()[1]-1) = 0.0;
        for (atlas::idx_t jl = fldInView.shape()[1]-2; jl > -1; --jl) {
          fldInView(jn, jl+1) += 0.5 * fldOutView(jn, jl);
          fldInView(jn, jl) += 0.5 * fldOutView(jn, jl);
          fldOutView(jn, jl) = 0.0;
        }
      }
    }
  }
}


}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VertLocInterp>
  makerVertLocInterp_("mo vertical interpolation for localization");

// -----------------------------------------------------------------------------
// Note that this assumes that the variables in activevars
//      exist in the variable level mappings of the inner and outer
//      oops Variables objects.
VertLocInterp::VertLocInterp(const oops::GeometryData & outerGeometryData,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const oops::FieldSet3D & xb,
                             const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    params_(params),
    outerGeometryData_(outerGeometryData),
    outerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    innerVars_(createInnerVars(params.innerVerticalLevels, activeVars_, outerVars))
{
  oops::Log::trace() << classname() << "::VertLocInterp starting" << std::endl;
  oops::Log::trace() << classname() << "::VertLocInterp done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLocInterp::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  // Takes active fields in fset that are on
  // half-level stagger and puts them on the correct vertical stagger
  // consistent with their variable names.

  // Create fieldset on functionspace
  atlas::FieldSet fsetOut;
  for (const auto & var : activeVars_) {
    // get levels from outer vars (active vars ambiguous)
    atlas::Field field =
      outerGeometryData_.functionSpace().createField<double>(
          atlas::option::name(var.name()) |
          atlas::option::levels(outerVars_[var.name()].getLevels()) |
          atlas::option::halo(1));
    atlas::array::make_view<double, 2>(field).assign(0.0);
    field.set_dirty(false);
    fsetOut.add(field);
  }

  // Simple vertical interpolation
  params_.reproduceBugStaggerDefn ?
    mo_buggy_verticalInterpolationFromHalfLevels(activeVars_, fset.fieldSet(), fsetOut)
    : verticalInterpolationFromFullLevels(activeVars_, fset.fieldSet(), fsetOut);

  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiply done"
                     << std::endl;
}

// -----------------------------------------------------------------------------

void VertLocInterp::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname()
                     << "::multiplyAD starting" << std::endl;

  // On input: fset on vertical staggers consistent with
  //           variable names

  // Create fields that are all on half-levels
  atlas::FieldSet fsetOut;
  for (const auto & var : activeVars_) {
    // get levels from inner vars (active vars ambiguous)
    atlas::Field field =
      outerGeometryData_.functionSpace().createField<double>
        (atlas::option::name(var.name()) |
         atlas::option::levels(innerVars_[var.name()].getLevels()) |
         atlas::option::halo(1));
    atlas::array::make_view<double, 2>(field).assign(0.0);
    field.set_dirty(false);
    fsetOut.add(field);
  }

  // Adjoint of simple vertical interpolation scheme
  params_.reproduceBugStaggerDefn ?
    mo_buggy_verticalInterpolationFromHalfLevelsAD(activeVars_, fsetOut, fset.fieldSet())
    : verticalInterpolationFromFullLevelsAD(activeVars_, fsetOut, fset.fieldSet());

  // keep passive vars
  for (auto & fld : fset) {
    if (!activeVars_.has(fld.name())) {
      fsetOut.add(fld);
    }
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname()
                     << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLocInterp::print(std::ostream & os) const {
  os << classname();
}

}  // namespace vader
}  // namespace saber
