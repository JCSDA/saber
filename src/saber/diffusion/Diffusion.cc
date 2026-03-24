/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/diffusion/Diffusion.h"

namespace saber {

// --------------------------------------------------------------------------------------

static SaberCentralBlockMaker<Diffusion> makerSaberDiffusion_("diffusion");

// --------------------------------------------------------------------------------------

Diffusion::Diffusion(
    const oops::GeometryData & geometryData,
    const oops::Variables & centralVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberCentralBlockBase(params, xb.validTime(), geometryData, centralVars),
    params_(params)
{
  // Compute total number of levels
  size_t nlevs = 0;
  for (const auto & var : params.activeVars.value().get_value_or(centralVars)) {
    nlevs += var.getLevels();
  }
  // Compute control vector size
  ctlVecSize_ = nlevs*geometryData.functionSpace().size();
}

// --------------------------------------------------------------------------------------

}  // namespace saber
