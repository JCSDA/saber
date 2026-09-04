/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierCovariance.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<BifourierCovariance> makerBifourierCovariance_("BifourierCovariance");

// -----------------------------------------------------------------------------

BifourierCovariance::BifourierCovariance(const oops::GeometryData & geometryData,
                                         const oops::Variables & centralVars,
                                         const eckit::Configuration & covarConf,
                                         const Parameters_ & params,
                                         const oops::FieldSet3D & xb,
                                         const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime(), geometryData, centralVars)
{
  oops::Log::trace() << classname() << "::BifourierCovariance starting" << std::endl;

  // Setup covariance implementation
  covar_ = std::make_unique<BifourierCovarianceImpl>(geometryData, centralVars, covarConf, params,
    xb, fg);

  oops::Log::trace() << classname() << "::BifourierCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierCovariance::~BifourierCovariance() {
  oops::Log::trace() << classname() << "::~BifourierCovariance starting" << std::endl;

  // Reset covariance implementation
  covar_.reset();

  oops::Log::trace() << classname() << "::~BifourierCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

