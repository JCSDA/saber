/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralAnalyticalCorrelation.h"

#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<SpectralAnalyticalCorrelation> makerSpectralAnalyticalCorrelation_(
        "spectral analytical filter");

// -----------------------------------------------------------------------------

SpectralAnalyticalCorrelation::SpectralAnalyticalCorrelation(
                                                   const oops::GeometryData & geometryData,
                                                   const oops::Variables & centralVars,
                                                   const eckit::Configuration & covarConf,
                                                   const Parameters_ & params,
                                                   const oops::FieldSet3D & xb,
                                                   const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime(), geometryData, centralVars)
{
  oops::Log::trace() << classname() << "::SpectralAnalyticalCorrelation starting " << std::endl;

  // Setup implementation
  impl_ = std::make_unique<SpectralAnalyticalImpl>(geometryData, centralVars, covarConf, params,
    xb.validTime(), true);

  oops::Log::trace() << classname() << "::SpectralAnalyticalCorrelation done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalCorrelation::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Adjoint multiplication
  impl_->multiplyAD(fset);

  // Forward multiplication
  impl_->multiply(fset);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalCorrelation::multiplySqrt(const atlas::Field & cv,
                                                 oops::FieldSet3D & fset,
                                                 const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Copy from control vector to fieldset
  cvToFset(cv, fset, offset, centralVars());

  // Forward multiplication
  impl_->multiply(fset);

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalCorrelation::multiplySqrtAD(const oops::FieldSet3D & fset,
                                                   atlas::Field & cv,
                                                   const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;

  // Copy
  oops::FieldSet3D fsetCopy(fset);

  // Adjoint multiplication
  impl_->multiplyAD(fsetCopy);

  // Copy from fieldset to control vector
  fsetToCv(fsetCopy, cv, offset, centralVars());

  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
