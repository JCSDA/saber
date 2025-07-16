/*
 * (C) Crown Copyright 2023-2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/GpToHpm1.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/util/Config.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/for_each.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

namespace {

using atlas::array::make_view;
using atlas::idx_t;

oops::Variables removeOuterOnlyVar(const oops::Variables & vars) {
  oops::Variables innerVars(vars);
  innerVars -= innerVars["hydrostatic_pressure_levels_minus_one"];
  return innerVars;
}

// ------------------------------------------------------------------------------------------------
// Note this is a copy of evalHydrostaticPressureTL
void eval_hydrostatic_pressure_levels_minus_one_tl(atlas::FieldSet & incFlds,
                                                   const atlas::FieldSet & augStateFlds) {
  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_tl()] starting ..." << std::endl;

  const idx_t levels = incFlds["geostrophic_pressure_levels_minus_one"].shape(1);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  // Loop over columns because fields have different numbers of levels
  util::for_each_column(
    [=](ConstView uPIncViewCol, View hPIncViewCol) {
      for (idx_t jl = 0; jl < levels; ++jl) {
        hPIncViewCol(jl) = uPIncViewCol(jl);
      }
    },
    incFlds["unbalanced_pressure_levels_minus_one"],
    incFlds["hydrostatic_pressure_levels_minus_one"]);

  if (augStateFlds.has("interpolation_weights")) {
    // Bins Vertical regression matrix stored in one field
    // B = (vertical regression matrix bin_0)
    //     (vertical regression matrix bin_1)
    //     (          ...                   )
    //     (vertical regression matrix bin_m)
    // Since each matrix is square we can easily infer the bin index from the row index
    //
    // First index of vertRegView is bin_index * number of levels + level index,
    // the second is number of levels associated with matrix column.
    //
    // First index of interpolation_weights is horizontal index, the second is bin index here

    const auto vertRegView = make_view<const double, 2>(
      augStateFlds["vertical_regression_matrices"]);

    const idx_t nBins = augStateFlds["interpolation_weights"].shape(1);

    using View = atlas::array::LocalView<double, 1>;
    using ConstView = atlas::array::LocalView<const double, 1>;

    util::for_each_column(
      [=](ConstView interpWeightViewCol, ConstView gPIncViewCol, View hPIncViewCol) {
        for (idx_t jl = 0; jl < levels; ++jl) {
          for (idx_t b = 0; b < nBins; ++b) {
            for (idx_t jl2 = 0; jl2 < levels; ++jl2) {
               hPIncViewCol(jl) += interpWeightViewCol(b) *
                                   vertRegView(b * levels + jl, jl2) *
                                   gPIncViewCol(jl2);
            }
          }
        }
      },
      augStateFlds["interpolation_weights"],
      incFlds["geostrophic_pressure_levels_minus_one"],
      incFlds["hydrostatic_pressure_levels_minus_one"]);
  }

  incFlds["hydrostatic_pressure_levels_minus_one"].set_dirty();

  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_tl()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
// Note this is a copy of evalHydrostaticPressureAD
void eval_hydrostatic_pressure_levels_minus_one_ad(atlas::FieldSet & hatFlds,
                                                   const atlas::FieldSet & augStateFlds) {
  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_ad()] starting ..." << std::endl;

  const idx_t levels = hatFlds["geostrophic_pressure_levels_minus_one"].shape(1);

  if (augStateFlds.has("interpolation_weights")) {
    // Bins Vertical regression matrix stored in one field
    // B = (vertical regression matrix bin_0)
    //     (vertical regression matrix bin_1)
    //     (          ...                   )
    //     (vertical regression matrix bin_m)
    // Since each matrix is square we can easily infer the bin index from the row index
    //
    // First index of vertRegView is bin_index * number of levels + level index,
    //     the second is level index

    const auto vertRegView = make_view<const double, 2>(
      augStateFlds["vertical_regression_matrices"]);

    const idx_t nBins = augStateFlds["interpolation_weights"].shape(1);

    using View = atlas::array::LocalView<double, 1>;
    using ConstView = atlas::array::LocalView<const double, 1>;

    util::for_each_column(
      [=](ConstView interpWeightViewCol,
          ConstView hPHatViewCol,
          View gPHatViewCol) {
        for (idx_t jl = levels - 1; jl >= 0; --jl) {
          for (idx_t b = nBins -1; b >= 0; --b) {
            for (idx_t jl2 = levels - 1; jl2 >= 0; --jl2) {
               gPHatViewCol(jl2) += interpWeightViewCol(b) *
                                    vertRegView(b * levels + jl, jl2) *
                                    hPHatViewCol(jl);
            }
          }
        }
      },
      augStateFlds["interpolation_weights"],
      hatFlds["hydrostatic_pressure_levels_minus_one"],
      hatFlds["geostrophic_pressure_levels_minus_one"]);
  }

  using View = atlas::array::LocalView<double, 1>;

  // Loop over columns because fields have different numbers of levels
  util::for_each_column(
    [=](View uPHatViewCol, View hPHatViewCol) {
      for (idx_t jl = levels - 1; jl >= 0; --jl) {
        uPHatViewCol(jl) += hPHatViewCol(jl);
        hPHatViewCol(jl) = 0.0;
      }
    },
    hatFlds["unbalanced_pressure_levels_minus_one"],
    hatFlds["hydrostatic_pressure_levels_minus_one"]);

  hatFlds["geostrophic_pressure_levels_minus_one"].set_dirty();
  hatFlds["unbalanced_pressure_levels_minus_one"].set_dirty();
  hatFlds["hydrostatic_pressure_levels_minus_one"].set_dirty();

  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_ad()] ... done" << std::endl;
}

// ------------------------------------------------------------------------------------------------
void eval_hydrostatic_pressure_levels_minus_one_tl_inv(atlas::FieldSet & incFlds,
                                                       const atlas::FieldSet & augStateFlds) {
  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_tl_inv()] starting ..." << std::endl;

  const idx_t levels = incFlds["unbalanced_pressure_levels_minus_one"].shape(1);

  using View = atlas::array::LocalView<double, 1>;
  using ConstView = atlas::array::LocalView<const double, 1>;

  // Loop over columns because fields have different numbers of levels
  util::for_each_column(
    [=](ConstView hPIncViewCol, View uPIncViewCol) {
      for (idx_t jl = 0; jl < levels; ++jl) {
        uPIncViewCol(jl) = hPIncViewCol(jl);
      }
    },
    incFlds["hydrostatic_pressure_levels_minus_one"],
    incFlds["unbalanced_pressure_levels_minus_one"]);

  if (augStateFlds.has("interpolation_weights")) {
    // Bins Vertical regression matrix stored in one field
    // B = (vertical regression matrix bin_0)
    //     (vertical regression matrix bin_1)
    //     (          ...                   )
    //     (vertical regression matrix bin_m)
    // Since each matrix is square we can easily infer the bin index from the row index
    //
    // First index of vertRegView is bin_index * number of levels + level index,
    // the second is number of levels associated with matrix column.
    //
    // First index of interpolation_weights is horizontal index, the second is bin index here

    const auto vertRegView = make_view<const double, 2>(
      augStateFlds["vertical_regression_matrices"]);

    const idx_t nBins = augStateFlds["interpolation_weights"].shape(1);

    using View = atlas::array::LocalView<double, 1>;
    using ConstView = atlas::array::LocalView<const double, 1>;

    util::for_each_column(
      [=](ConstView interpWeightViewCol, ConstView gPIncViewCol, View uPIncViewCol) {
        for (idx_t jl = 0; jl < levels; ++jl) {
          for (idx_t b = 0; b < nBins; ++b) {
            for (idx_t jl2 = 0; jl2 < levels; ++jl2) {
               uPIncViewCol(jl) -= interpWeightViewCol(b) *
                                   vertRegView(b * levels + jl, jl2) *
                                   gPIncViewCol(jl2);
            }
          }
        }
      },
      augStateFlds["interpolation_weights"],
      incFlds["geostrophic_pressure_levels_minus_one"],
      incFlds["unbalanced_pressure_levels_minus_one"]);
  }

  incFlds["unbalanced_pressure_levels_minus_one"].set_dirty();

  oops::Log::trace()
    << "[eval_hydrostatic_pressure_levels_minus_one_tl_inv()] ... done" << std::endl;
}


}  // namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GpToHpm1>
  makerGpToHpm1_("mo_hydrostatic_pressure_levels_minus_one_from_geostrophic_pressure");

// -----------------------------------------------------------------------------

GpToHpm1::GpToHpm1(const oops::GeometryData & outerGeometryData,
                   const oops::Variables & outerVars,
                   const eckit::Configuration & covarConf,
                   const Parameters_ & params,
                   const oops::FieldSet3D & xb,
                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(removeOuterOnlyVar(getUnionOfInnerActiveAndOuterVars(params, outerVars))),
    activeOuterVars_(params.activeOuterVars(outerVars)),
    innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
    params_(params),
    covFieldSet_(),
    augmentedStateFieldSet_()
{
  oops::Log::trace() << classname() << "::GpToHpm1 starting" << std::endl;
  const oops::Variables stateVariables = params.mandatoryStateVars();
  augmentedStateFieldSet_.clear();
  for (const auto & s : stateVariables) {
    augmentedStateFieldSet_.add(xb.fieldSet()[s.name()]);
  }

  oops::Log::trace() << classname() << "::GpToHpm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

GpToHpm1::~GpToHpm1() {
  oops::Log::trace() << classname() << "::~GpToHpm1 starting" << std::endl;
  util::Timer timer(classname(), "~GpToHpm1");
  oops::Log::trace() << classname() << "::~GpToHpm1 done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  // Allocate output fields if they are not already present, e.g when randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  // Populate output fields.
  eval_hydrostatic_pressure_levels_minus_one_tl(fset.fieldSet(),
                                                augmentedStateFieldSet_);

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  eval_hydrostatic_pressure_levels_minus_one_ad(fset.fieldSet(),
                                                augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::multiplyAD done"  << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  if (!fset.has("geostrophic_pressure_levels_minus_one")) {
    oops::Log::error() << "The inverse of block "
          << classname()
          << " is not correctly defined if geostrophic_pressure_levels_minus_one"
          << " is not provided as an input." << std::endl;
    throw eckit::UserError("Please only use leftInverseMultiply of this block "
                           "within the mo_hydrostatic_pressure_levels_minus_one "
                           "block.", Here());
  }
  //   Allocate inner-only variables except air temperature
  oops::Variables innerOnlyVarsForInversion(innerOnlyVars_);
  innerOnlyVarsForInversion -=
    innerOnlyVarsForInversion["geostrophic_pressure_levels_minus_one"];
  checkFieldsAreNotAllocated(fset, innerOnlyVarsForInversion);
  allocateMissingFields(fset, innerOnlyVarsForInversion, innerOnlyVarsForInversion,
                        innerGeometryData_.functionSpace());

  // Retrieve unbalanced pressure from hydrostatic pressure and geostrophic pressure.
  eval_hydrostatic_pressure_levels_minus_one_tl_inv(fset.fieldSet(),
                                                    augmentedStateFieldSet_);
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

void GpToHpm1::read() {
  oops::Log::trace() << classname() << "::read start " << params_ <<  std::endl;
  const auto & readParams = params_.readParams.value();
  if (readParams != boost::none) {
    eckit::LocalConfiguration lconf;
    readParams.value().serialize(lconf);
    const eckit::Configuration & conf = lconf;
    // Covariance FieldSet
    covFieldSet_ = createGpRegressionStats(innerGeometryData_.functionSpace(),
                                           innerVars_,
                                           conf);
    // also copy variables from covariance fieldset if required
    if (covFieldSet_.has("interpolation_weights")) {
      augmentedStateFieldSet_.add(covFieldSet_["vertical_regression_matrices"]);
      augmentedStateFieldSet_.add(covFieldSet_["interpolation_weights"]);
    }
  }
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

void GpToHpm1::directCalibration(const oops::FieldSets & fset) {
  oops::Log::info() << classname() << "::directCalibration start" << std::endl;
  const auto & calibrationReadParams = params_.calibrationReadParams.value();
  if (calibrationReadParams != boost::none) {
    eckit::LocalConfiguration lconf;
    calibrationReadParams.value().serialize(lconf);
    const eckit::Configuration & conf = lconf;
    // Covariance FieldSet
    covFieldSet_ = createGpRegressionStats(innerGeometryData_.functionSpace(),
                                           innerVars_,
                                           conf);

    // also copy variables from covariance fieldset if required
    if (covFieldSet_.has("interpolation_weights")) {
      augmentedStateFieldSet_.add(covFieldSet_["vertical_regression_matrices"]);
      augmentedStateFieldSet_.add(covFieldSet_["interpolation_weights"]);
    }
  }
  oops::Log::info() << classname() << "::directCalibration end" << std::endl;
}

void GpToHpm1::write() const {
  oops::Log::trace() << classname() << "::write start" << std::endl;
  // write regression matrix to file.
  oops::Log::trace() << classname() << "::write end" << std::endl;
}

// -----------------------------------------------------------------------------

void GpToHpm1::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
