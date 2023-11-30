/*
 * (C) Crown Copyright 2023 Met Office
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
    if (!function.has("horizontal daley length")) throw eckit::BadParameter(
                    "horizontal daley length must be specified for Gaussian function shapes");

    const double loc_rh = function.getDouble("horizontal daley length");
    gaussianShape(spectralFilter, loc_rh);
  } else {
    throw eckit::BadParameter("function shape " + functionShape + " not implemented yet.",
                              Here());
  }

  // 3) Normalize as a localization function if required
  // ---------------------------------------------------
  if (params.normalize.value()) {
    // Compute total variance in spectral space before normalization.
    double totalVarianceSpectral = 0;
    const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
    for (int jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
      const int m = zonal_wavenumbers(jm);
      for (int n = m; n <= truncation; ++n) {
        if (m == 0) {
          totalVarianceSpectral += spectralFilter[n];
        } else {
          // both real and imaginary components
          totalVarianceSpectral += 2 * spectralFilter[n];
        }
      }
    }
    geometryData.comm().allReduceInPlace(totalVarianceSpectral, eckit::mpi::sum());

    // Normalize total variance to be 1
    const double normalization = 1.0 / totalVarianceSpectral;
    std::transform(spectralFilter.begin(), spectralFilter.end(),
                   spectralFilter.begin(),
                   [& normalization](auto & elem){return elem * normalization;});
  }

  // 4) Take square root
  // -------------------
  std::transform(spectralFilter.begin(), spectralFilter.end(),
                 spectralFilter.begin(), [](auto & e){return std::sqrt(e);});

  oops::Log::trace() << "saber::spectralb::createSpectralFilter done"
                     << std::endl;
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
    activeVars_(getActiveVars(params, outerVars)),
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

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const int truncation = specFunctionSpace_.truncation();
  const int max_zonal_wavenumber = zonal_wavenumbers.size();

  for (const auto & var : activeVars_.variables()) {
    auto view = atlas::array::make_view<double, 2>(fieldSet[var]);
    const atlas::idx_t levels = fieldSet[var].levels();

    // Element-wise multiplication
    atlas::idx_t jnode = 0;
    for (int jm = 0; jm < max_zonal_wavenumber; ++jm) {
      int m = zonal_wavenumbers(jm);
      for (int n = m; n <= truncation; ++n) {
        // Real component
        for (atlas::idx_t jlev = 0; jlev < levels; jlev++) {
          view(jnode, jlev) *= spectralFilter_[n];
        }
        jnode++;
        // Imaginary component
        for (atlas::idx_t jlev = 0; jlev < levels; jlev++) {
          view(jnode, jlev) *= spectralFilter_[n];
        }
        jnode++;
      }
    }
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

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const int truncation = specFunctionSpace_.truncation();
  const int max_zonal_wavenumber = zonal_wavenumbers.size();

  std::vector<double> filterInverse(spectralFilter_.size());
  std::transform(spectralFilter_.begin(), spectralFilter_.end(),
                 filterInverse.begin(), [&](auto & x){return x != 0.0 ? 1.0 / x : 0.0;});

  for (const auto & var : activeVars_.variables()) {
    auto view = atlas::array::make_view<double, 2>(fieldSet[var]);
    const atlas::idx_t levels = fieldSet[var].levels();

    // Element-wise multiplication
    atlas::idx_t jnode = 0;
    for (int jm = 0; jm < max_zonal_wavenumber; ++jm) {
      int m = zonal_wavenumbers(jm);
      for (int n = m; n <= truncation; ++n) {
        // Real component
        for (atlas::idx_t jlev = 0; jlev < levels; jlev++) {
          view(jnode, jlev) *= filterInverse[n];
        }
        jnode++;
        // Imaginary component
        for (atlas::idx_t jlev = 0; jlev < levels; jlev++) {
          view(jnode, jlev) *= filterInverse[n];
        }
        jnode++;
      }
    }
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
