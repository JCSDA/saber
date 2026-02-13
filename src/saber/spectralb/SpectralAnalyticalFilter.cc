/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

#include "saber/spectralb/SpectralAnalyticalFilter.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Earth.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Operation.h"

#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

#include "saber/oops/Utilities.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SpectralAnalyticalFilter> makerSpectralAnalyticalFilter_(
        "spectral analytical filter");

// -----------------------------------------------------------------------------

namespace {

void gaussianShape(std::vector<double> & wavenumbers,
                   const double loc_rh) {
  // Transform wavenumbers into Gaussian
  oops::Log::info() << "Info     : Horizontal localization length is "
                    << loc_rh / 1000.0 << " km." << std::endl;
  const double standard_wavenumber = atlas::util::Earth::radius() / loc_rh;
  const double factor = 1.0 / (2 * std::pow(standard_wavenumber, 2));
  std::transform(wavenumbers.begin(), wavenumbers.end(),
                 wavenumbers.begin(),
                 [& factor](auto & n){return std::exp(- factor * std::pow(n, 2));});
}

void triangularShape(std::vector<double> & wavenumbers,
                     const double & bandmin, const double & bandpeak,
                     const double & bandmax) {
  for (int w = 0; w < static_cast<int>(wavenumbers.size()); ++w) {
    double left = -(bandmin - w) / (bandpeak - bandmin);
    double right = -(bandmax - w) / (bandpeak - bandmax);
    wavenumbers[w] = std::max(0.0, std::min(left, right));
  }
}

void rectangleRightAngleTriangleShape(std::vector<double> & wavenumbers,
    const double & bandmin, const double & bandpeak, const double & bandmax) {
  for (int w = 0; w < static_cast<int>(wavenumbers.size()); ++w) {
    double left = 1.0;
    double right = -(bandmax - w) / (bandpeak - bandmax);
    wavenumbers[w] = std::max(0.0, std::min(left, right));
  }
}

void rightAngleTriangleRectangleShape(std::vector<double> & wavenumbers,
    const double & bandmin, const double & bandpeak, const double & bandmax) {
  for (int w = 0; w < static_cast<int>(wavenumbers.size()); ++w) {
    double left = -(bandmin - w) / (bandpeak - bandmin);
    double right = 1.0;
    wavenumbers[w] = std::max(0.0, std::min(left, right));
  }
}

void rectangleShape(std::vector<double> & wavenumbers,
    const double & bandmin, const double & amplitude, const double & bandmax) {
  for (int w = 0; w < static_cast<int>(wavenumbers.size()); ++w) {
    wavenumbers[w] = (w >= bandmin && w <= bandmax ? amplitude : 0.0);
  }
}

// -----------------------------------------------------------------------------

auto createSpectralFilter(const oops::GeometryData & geometryData,
                          const SpectralAnalyticalFilterParameters & params) {
  oops::Log::trace() << "saber::spectralb::createSpectralFilter starting"
                     << std::endl;

  // 1) Initialize filter as a vector of total wavenumbers
  // -----------------------------------------------------
  const atlas::functionspace::Spectral specFunctionSpace(geometryData.functionSpace());
  const int truncation = specFunctionSpace.truncation();
  std::vector<double> spectralFilter(truncation + 1);
  std::iota(spectralFilter.begin(), spectralFilter.end(), 0);

  // 2) Modify into required filter
  // ------------------------------
  const auto & function = params.function.value();
  const std::string functionShape(function.getString("shape", "gaussian"));
  if (functionShape.compare("gaussian") == 0) {
    if (!function.has("horizontal daley length")) {
      throw eckit::BadParameter(
      "horizontal daley length must be specified for Gaussian function shapes");
    }
    const double loc_rh = function.getDouble("horizontal daley length");
    gaussianShape(spectralFilter, loc_rh);
  } else if (functionShape.compare("waveband filter") == 0) {
    if (!function.has("waveband min")) throw eckit::BadParameter(
      "minimum wavenumber must be specified for waveband function shapes");
    if (!function.has("waveband max")) throw eckit::BadParameter(
      "waveband max wavenumber must be specified for waveband function shapes");
    if (!function.has("waveband peak")) throw eckit::BadParameter(
      "waveband peak wavenumber must be specified for waveband function shapes");
    const int bandmin = function.getInt("waveband min");
    const int bandmax = function.getInt("waveband max");
    const int bandpeak = function.getInt("waveband peak");

    (bandmin == 0 ?
     rectangleRightAngleTriangleShape(spectralFilter, 0.0,
                                      static_cast<double>(bandpeak),
                                      static_cast<double>(bandmax)) :
       (bandmax == truncation ?
        rightAngleTriangleRectangleShape(spectralFilter,
                                         static_cast<double>(bandmin),
                                         static_cast<double>(bandpeak),
                                         static_cast<double>(truncation)) :
        triangularShape(spectralFilter,
                        static_cast<double>(bandmin),
                        static_cast<double>(bandpeak),
                        static_cast<double>(bandmax))));
  } else if (functionShape.compare("boxcar") == 0) {
    if (!function.has("waveband min")) throw eckit::BadParameter(
      "minimum wavenumber must be specified for boxcar functions");
    if (!function.has("waveband max")) throw eckit::BadParameter(
      "maximum wavenumber must be specified for boxcar functions");
    if (!function.has("waveband amplitude")) throw eckit::BadParameter(
      "waveband amplitude must be specified for boxcar functions");
    const int bandmin = function.getInt("waveband min");
    const int bandmax = function.getInt("waveband max");
    const int amplitude = function.getInt("waveband amplitude");
    rectangleShape(spectralFilter, bandmin, amplitude, bandmax);
  } else {
    throw eckit::BadParameter("function shape " + functionShape + " not implemented yet.",
                              Here());
  }

  // 3) Normalize as a localization function if required or
  //    preserve variance of increments.
  // ---------------------------------------------------
  if ( params.normalizeFilterVariance && params.preservingVariance ) {
    throw eckit::BadParameter(
    "normalize filter variance option incompatibile with preserving variance option");
  }

  if ( params.normalizeFilterVariance ) {
    // Compute total variance in spectral space before normalization.
    std::vector<double> spectralFilterPrefixSum(spectralFilter.size(), 0.0);
    std::partial_sum(spectralFilter.begin(), spectralFilter.end(), spectralFilterPrefixSum.begin());

    double totalVarianceSpectral = 0;
    const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
    for (int jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
      const int m = zonal_wavenumbers(jm);
      if (m == 0) {
        totalVarianceSpectral += spectralFilterPrefixSum[truncation];
      } else if (m <= truncation) {
        // both real and imaginary components
        totalVarianceSpectral += 2 * (spectralFilterPrefixSum[truncation]
                                  - spectralFilterPrefixSum[m-1]);
      }
    }
    geometryData.comm().allReduceInPlace(totalVarianceSpectral, eckit::mpi::sum());

    // Normalize total variance to be 1
    const double normalization = 1.0 / totalVarianceSpectral;
    std::transform(spectralFilter.begin(), spectralFilter.end(),
                   spectralFilter.begin(),
                   [& normalization](auto & elem){return elem * normalization;});
  }

  if ( params.preservingVariance ) {
    std::transform(spectralFilter.begin(), spectralFilter.end(),
                   spectralFilter.begin(), [](auto & e){return std::sqrt(e);});
  }

  if ( params.complementFilter ) {
    std::transform(spectralFilter.begin(), spectralFilter.end(),
                   spectralFilter.begin(), [](auto & e){return 1.0 - e;});
  }

  // 4) Take square root (as this is an outer block)
  // -----------------------------------------------
  std::transform(spectralFilter.begin(), spectralFilter.end(),
                 spectralFilter.begin(), [](auto & e){return std::sqrt(e);});

  return spectralFilter;
}

}  // namespace

// -----------------------------------------------------------------------------
SpectralAnalyticalFilter::SpectralAnalyticalFilter(const oops::GeometryData & geometryData,
                                                   const oops::Variables & outerVars,
                                                   const eckit::Configuration & covarConf,
                                                   const Parameters_ & params,
                                                   const oops::FieldSet3D & xb,
                                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()), params_(params),
    activeVars_(params.getActiveVars(outerVars)),
    innerGeometryData_(geometryData),
    innerVars_(outerVars),
    specFunctionSpace_(geometryData.functionSpace()),
    spectralFilter_(createSpectralFilter(geometryData, params))
{
  oops::Log::trace() << classname() << "::SpectralAnalyticalFilter starting " << std::endl;
  oops::Log::trace() << classname() << "::SpectralAnalyticalFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalFilter::multiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  for (const auto & var : activeVars_) {
    auto view = atlas::array::make_view<double, 2>(fieldSet[var.name()]);

    // Element-wise multiplication
    util::for_each_index(
      util::make_index_space_spectral_1d(fieldSet[var.name()]),
      [&](atlas::idx_t i_part, atlas::idx_t n, atlas::idx_t m) mutable {
        for (int l = 0; l < view.shape(1); ++l) {
          view(i_part, l) *= spectralFilter_[n];
        }
      });
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalFilter::multiplyAD(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // The block is self-adjoint:
  multiply(fieldSet);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalFilter::leftInverseMultiply(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  std::vector<double> filterInverse(spectralFilter_.size());
  std::transform(spectralFilter_.begin(), spectralFilter_.end(),
                 filterInverse.begin(), [&](auto & x){return x != 0.0 ? 1.0 / x : 0.0;});

  for (const auto & var : activeVars_) {
    auto view = atlas::array::make_view<double, 2>(fieldSet[var.name()]);

    // Element-wise multiplication
    util::for_each_index(
      util::make_index_space_spectral_1d(fieldSet[var.name()]),
      [&](atlas::idx_t i_part, atlas::idx_t n, atlas::idx_t m) mutable {
        for (int l = 0; l < view.shape(1); ++l) {
          view(i_part, l) *= filterInverse[n];
        }
      });
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D SpectralAnalyticalFilter::generateInnerFieldSet(
  const oops::GeometryData & innerGeometryData,
  const oops::Variables & innerVars) const {
  oops::FieldSet3D fset(this->validTime(), innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D SpectralAnalyticalFilter::generateOuterFieldSet(
  const oops::GeometryData & outerGeometryData,
  const oops::Variables & outerVars) const {
  oops::FieldSet3D fset(this->validTime(), outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void SpectralAnalyticalFilter::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
