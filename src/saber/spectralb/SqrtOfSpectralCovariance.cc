/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SqrtOfSpectralCovariance.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"
#include "saber/spectralb/CovarianceStatisticsUtils.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SqrtOfSpectralCovariance>
    makerSqrtOfSpectralCovariance_("square root of spectral covariance");

// -----------------------------------------------------------------------------
SqrtOfSpectralCovariance::SqrtOfSpectralCovariance(
    const oops::GeometryData & outerGeometryData,
    const oops::Variables & outerVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    params_(params),
    activeVars_(getActiveVars(params, outerVars)),
    outerVars_(outerVars),
    specFunctionSpace_(outerGeometryData.functionSpace()),
    innerGeometryData_(outerGeometryData)
{
  oops::Log::trace() << classname() << "::SqrtOfSpectralCovariance starting " << std::endl;

  oops::Log::trace() << classname() << "::SqrtOfSpectralCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  specutils::spectralVerticalConvolutionSqrt(activeVars_,
                                             specFunctionSpace_,
                                             spectralUMatrices_,
                                             fieldSet.fieldSet());

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyUMatrixAD starting" << std::endl;

  specutils::spectralVerticalConvolutionSqrtAD(activeVars_,
                                               specFunctionSpace_,
                                               spectralUMatrices_,
                                               fieldSet.fieldSet());

  oops::Log::trace() << classname() << "::multiplyUMatrixAD done" << std::endl;
}



// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Note this is reading in the full spectral bins that exist in the cov file
  spectralUMatrices_.clear();
  const spectralbReadParameters readP = *params_.readParams.value();
  const std::vector<std::size_t> nSpectralBinsFull =
    specutils::getNSpectralBinsFull(readP, activeVars_);
  spectralUMatrices_ =
    specutils::createUMatrices(activeVars_, nSpectralBinsFull, readP);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
