/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "saber/spectralb/SqrtOfSpectralCovariance.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"

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
    const atlas::FieldSet & xb,
    const atlas::FieldSet & fg,
    const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    params_(params),
    activeVars_(getActiveVars(params, outerVars)),
    outerVars_(outerVars),
    cs_(),
    specFunctionSpace_(outerGeometryData.functionSpace()),
    innerGeometryData_(outerGeometryData)
{
  oops::Log::trace() << classname() << "::SqrtOfSpectralCovariance starting " << std::endl;

  if (params.doCalibration()) {
    // TODO(Marek)
  } else {
    variance_opt_ = params.readParams.value()->varianceOpt;
  }

  oops::Log::trace() << classname() << "::SqrtOfSpectralCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  using atlas::array::make_view;
  using atlas::idx_t;

  const idx_t N = specFunctionSpace_.truncation();

  // 1D convolution with spectral vertical covariances for each total wavenumber n1
  const atlas::FieldSet & spectralUMatrices = (variance_opt_ ?
                                         cs_->getSpectralUMatrices() :
                                         cs_->getSpectralCorrelUMatrices());

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t levels(fieldSet[var].levels());
    auto UMatrixView = make_view<const double, 3>(spectralUMatrices[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);
    const int nSpectralBinsFull = spectralUMatrices[var].shape(0);

    idx_t i = 0;
    std::vector<double> col(levels), col2(levels);
    for (idx_t jm1 = 0; jm1 < nb_zonal_wavenumbers; ++jm1) {
      const idx_t m1 = zonal_wavenumbers(jm1);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
        // note that img stands for imaginary component are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[static_cast<std::size_t>(jl)] = spfView(i, jl);
          }
          for (idx_t r = 0; r < levels; ++r) {
            col2[static_cast<std::size_t>(r)] = 0;
            for (idx_t c = 0; c < levels; ++c) {
              col2[static_cast<std::size_t>(r)] += UMatrixView(n1, r, c) * col[c];
            }
          }
          const double norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[static_cast<std::size_t>(jl)] / norm;
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyUMatrixAD starting" << std::endl;

  using atlas::array::make_view;
  using atlas::idx_t;

  const idx_t N = specFunctionSpace_.truncation();

  // 1D convolution with spectral vertical covariances for each total wavenumber n1
  const atlas::FieldSet & spectralUMatrices = (variance_opt_ ?
                                         cs_->getSpectralUMatrices() :
                                         cs_->getSpectralCorrelUMatrices());

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t levels(fieldSet[var].levels());
    auto UMatrixView = make_view<const double, 3>(spectralUMatrices[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);
    const int nSpectralBinsFull = spectralUMatrices[var].shape(0);

    idx_t i = 0;
    std::vector<double> col(levels), col2(levels);
    for (idx_t jm1 = 0; jm1 < nb_zonal_wavenumbers; ++jm1) {
      const idx_t m1 = zonal_wavenumbers(jm1);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
        // note that img stands for imaginary component are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[static_cast<std::size_t>(jl)] = spfView(i, jl);
          }
          for (idx_t r = 0; r < levels; ++r) {
            col2[static_cast<std::size_t>(r)] = 0;
            for (idx_t c = 0; c < levels; ++c) {
              col2[static_cast<std::size_t>(r)] += UMatrixView(n1, c, r) * col[c];
            }
          }
          const double norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[static_cast<std::size_t>(jl)] / norm;
          }
        }
      }
    }
  }
  oops::Log::trace() << classname() << "::multiplyUMatrixAD done" << std::endl;
}



// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Initialize CovStat_ErrorCov
  cs_.reset(new CovStat_ErrorCov(activeVars_,
                                 *params_.readParams.value()));


  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void SqrtOfSpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
