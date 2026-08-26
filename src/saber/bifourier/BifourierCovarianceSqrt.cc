/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierCovarianceSqrt.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierCovarianceSqrt>
  makerBifourierCovarianceSqrt_("BifourierCovarianceSqrt");

// -----------------------------------------------------------------------------

BifourierCovarianceSqrt::BifourierCovarianceSqrt(const oops::GeometryData & outerGeometryData,
                                                 const oops::Variables & outerVars,
                                                 const eckit::Configuration & covarConf,
                                                 const Parameters_ & params,
                                                 const oops::FieldSet3D & xb,
                                                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars)
{
  oops::Log::trace() << classname() << "::BifourierCovarianceSqrt starting" << std::endl;

  // Setup covariance implementation
  covar_ = std::make_unique<BifourierCovarianceImpl>(outerGeometryData, outerVars, covarConf,
    params, xb, fg);

  oops::Log::trace() << classname() << "::BifourierCovarianceSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierCovarianceSqrt::~BifourierCovarianceSqrt() {
  oops::Log::trace() << classname() << "::~BifourierCovarianceSqrt starting" << std::endl;

  // Reset covariance implementation
  covar_.reset();

  oops::Log::trace() << classname() << "::~BifourierCovarianceSqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
