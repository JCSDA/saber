/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralAnalyticalFilter.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SpectralAnalyticalFilter> makerSpectralAnalyticalFilter_(
        "spectral analytical filter");

// -----------------------------------------------------------------------------

SpectralAnalyticalFilter::SpectralAnalyticalFilter(const oops::GeometryData & geometryData,
                                                   const oops::Variables & outerVars,
                                                   const eckit::Configuration & covarConf,
                                                   const Parameters_ & params,
                                                   const oops::FieldSet3D & xb,
                                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), geometryData, outerVars)
{
  oops::Log::trace() << classname() << "::SpectralAnalyticalFilter starting " << std::endl;

  // Setup implementation
  impl_ = std::make_unique<SpectralAnalyticalImpl>(geometryData, outerVars, covarConf, params,
    xb.validTime(), false);

  oops::Log::trace() << classname() << "::SpectralAnalyticalFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
