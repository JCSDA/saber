/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/diffusion/DiffusionFilter.h"

namespace saber {

// --------------------------------------------------------------------------------------

static SaberOuterBlockMaker<DiffusionFilter> makerSaberDiffusionFilter_("diffusion filter");

// --------------------------------------------------------------------------------------

DiffusionFilter::DiffusionFilter(
    const oops::GeometryData & outerGeometryData,
    const oops::Variables & outerVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    params_(params)
{
  // Compute total number of levels
  size_t nlevs = 0;
  for (const auto & var : params.activeVars.value().get_value_or(outerVars)) {
    nlevs += var.getLevels();
  }
}
// --------------------------------------------------------------------------------------

}  // namespace saber
