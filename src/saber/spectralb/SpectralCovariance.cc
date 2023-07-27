/*
 * (C) Crown Copyright 2022-2023 Met Office
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "saber/spectralb/SpectralCovariance.h"

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

static SaberCentralBlockMaker<SpectralCovariance> makerSpectralCovariance_("spectral covariance");

// -----------------------------------------------------------------------------
SpectralCovariance::SpectralCovariance(const oops::GeometryData & geometryData,
                                       const oops::Variables & centralVars,
                                       const eckit::Configuration & covarConf,
                                       const Parameters_ & params,
                                       const atlas::FieldSet & xb,
                                       const atlas::FieldSet & fg,
                                       const util::DateTime & validTimeOfXbFg,
                                       const size_t & timeRank)
  : SaberCentralBlockBase(params), params_(params),
    activeVars_(getActiveVars(params, centralVars)),
    cs_(),
    geometryData_(geometryData),
    specFunctionSpace_(geometryData_.functionSpace()),
    timeRank_(timeRank)
{
  oops::Log::trace() << classname() << "::SpectralCovariance starting " << std::endl;

  if (params.doCalibration()) {
    // TODO(Marek)
  } else {
    variance_opt_ = params.readParams.value()->varianceOpt;
  }

  oops::Log::trace() << classname() << "::SpectralCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::randomize(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  oops::Log::error() << "randomization with spectral covariance saber block"
                     << " is not supported. Instead please use 'ID' central block"
                     << " and 'square root of spectral covariance' outer block."
                     << std::endl;
  throw(eckit::FunctionalityNotSupported(
        "use ID and square root of spectral covariance instead.", Here()));
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  using atlas::array::make_view;
  using atlas::idx_t;

  const idx_t N = specFunctionSpace_.truncation();

  const atlas::FieldSet & spectralVerticalCovariances = (variance_opt_ ?
                                         cs_->getSpectralVerticalCovariances() :
                                         cs_->getSpectralVerticalCorrelations());
  std::vector<std::string> vertCovNames = spectralVerticalCovariances.field_names();

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t i = 0;
    idx_t levels(fieldSet[var].levels());
    auto vertCovView = make_view<const double, 3>(spectralVerticalCovariances[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);

    std::vector<double> col(levels), col2(levels);
    // For each total wavenumber n1, perform a 1D convolution with vertical covariances.
    for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const idx_t m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
        // Note that img stands for imaginary component and are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img) {
          // Pre-fill vertical column to be convolved.
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[static_cast<std::size_t>(jl)] = spfView(i, jl);
          }
          // The 2*n1+1 factor is there to equally distribute the covariance across
          // the spectral coefficients associated to this total wavenumber.
          const double norm = static_cast<double>((2 * n1 + 1) * vertCovView.shape(0));
          for (idx_t r = 0; r < levels; ++r) {
            col2[static_cast<std::size_t>(r)] = 0;
            for (idx_t c = 0; c < levels; ++c) {
              col2[static_cast<std::size_t>(r)] += vertCovView(n1, r, c) * col[c] / norm;
            }
          }
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[static_cast<std::size_t>(jl)];
          }
          ++i;
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Initialize CovStat_ErrorCov
  cs_.reset(new CovStat_ErrorCov(activeVars_,
                                 *params_.readParams.value()));

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
