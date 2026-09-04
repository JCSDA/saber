/*
 * (C) Crown Copyright 2023-2024 Met Office
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

#include "saber/bifourier/BifourierAnalyticalFilter.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Earth.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Operation.h"

#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

#include "saber/bifourier/BifourierUtilities.h"
#include "saber/oops/Utilities.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierAnalyticalFilter> makerBifourierAnalyticalFilter_(
        "BifourierAnalyticalFilter");

// -----------------------------------------------------------------------------

namespace {

auto createSpectralFilter(const BifourierTransformBase & trans,
                          const BifourierAnalyticalFilterParameters & params) {
  oops::Log::trace() << "saber::bifourier::createSpectralFilter starting" << std::endl;

  // Initialize spectral filter
  std::vector<double> spectralFilter(trans.ns(), 0.0);

  // Get waveband peak
  const double wpeak = params.wavebandPeak.value();

  if (params.wavebandMin.value() && params.wavebandMax.value()) {
    // Intermediate scales

    // Get waveband bounds
    const double wmin = *params.wavebandMin.value();
    const double wmax = *params.wavebandMax.value();

    for (size_t js = 0; js < trans.ns(); ++js) {
      // Total wavenumber
      const double kstar = trans.rkstar(trans.k(js), trans.l(js), trans.M(), trans.N(),
        trans.nwGlb());

      // Response function
      if ((kstar > wmin) && (kstar <= wpeak)) {
        spectralFilter[js] = std::pow(std::sin(0.5*M_PI*(kstar-wmin)/(wpeak-wmin)), 2);
      } else if ((kstar > wpeak) && (kstar < wmax)) {
        spectralFilter[js] = std::pow(std::sin(0.5*M_PI*(wmax-kstar)/(wmax-wpeak)), 2);
      }
    }
  } else if (params.wavebandMin.value()) {
    // Smallest scales

    // Get waveband bounds
    const double wmin = *params.wavebandMin.value();

    for (size_t js = 0; js < trans.ns(); ++js) {
      // Total wavenumber
      const double kstar = trans.rkstar(trans.k(js), trans.l(js), trans.M(), trans.N(),
        trans.nwGlb());

      // Response function
      if (kstar > wpeak) {
        spectralFilter[js] = 1.0;
      } else if (kstar > wmin) {
        spectralFilter[js] = std::pow(std::sin(0.5*M_PI*(kstar-wmin)/(wpeak-wmin)), 2);
      }
    }
  } else if (params.wavebandMax.value()) {
    // Largest scales

    // Get waveband bounds
    const double wmax = *params.wavebandMax.value();

    for (size_t js = 0; js < trans.ns(); ++js) {
      // Total wavenumber
      const double kstar = trans.rkstar(trans.k(js), trans.l(js), trans.M(), trans.N(),
        trans.nwGlb());

      // Response function
      if (kstar < wpeak) {
        spectralFilter[js] = 1.0;
      } else if (kstar <= wmax) {
        spectralFilter[js] = std::pow(std::sin(0.5*M_PI*(wmax-kstar)/(wmax-wpeak)), 2);
      }
    }
  } else {
    // Missing min or max waveband
    throw eckit::BadParameter("min or max waveband should be present", Here());
  }

  if (params.inverseMode.value()) {
    // Inverse filter
    std::transform(spectralFilter.begin(), spectralFilter.end(),
                   spectralFilter.begin(), [&](auto & x){return x != 0.0 ? 1.0 / x : 0.0;});
  }

  oops::Log::trace() << "saber::bifourier::createSpectralFilter done" << std::endl;
  return spectralFilter;
}

}  // namespace

// -----------------------------------------------------------------------------
BifourierAnalyticalFilter::BifourierAnalyticalFilter(const oops::GeometryData & outerGeometryData,
                                                     const oops::Variables & outerVars,
                                                     const eckit::Configuration & covarConf,
                                                     const Parameters_ & params,
                                                     const oops::FieldSet3D & xb,
                                                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    activeVars_(params.getActiveVars(outerVars)),
    params_(params),
    trans_(transStore_.retrieveTransform(outerGeometryData)),
    spectralFilter_(createSpectralFilter(*trans_, params))
{
  oops::Log::trace() << classname() << "::BifourierAnalyticalFilter starting " << std::endl;
  oops::Log::trace() << classname() << "::BifourierAnalyticalFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAnalyticalFilter::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Get field views
    auto view = getView2D(var, fset);

    // Apply filter
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        view(js, jz) *= spectralFilter_[js];
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAnalyticalFilter::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // The block is self-adjoint:
  multiply(fset);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAnalyticalFilter::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  std::vector<double> filterInverse(spectralFilter_.size());
  std::transform(spectralFilter_.begin(), spectralFilter_.end(),
                 filterInverse.begin(), [&](auto & x){return x != 0.0 ? 1.0 / x : 0.0;});

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Get field views
    auto view = getView2D(var, fset);

    // Apply filter
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        view(js, jz) *= filterInverse[js];
      }
    }
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D BifourierAnalyticalFilter::generateInnerFieldSet(
  const oops::GeometryData & innerGeometryData,
  const oops::Variables & innerVars) const {
  oops::FieldSet3D fset(validTime_, innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));
  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D BifourierAnalyticalFilter::generateOuterFieldSet(
  const oops::GeometryData & outerGeometryData,
  const oops::Variables & outerVars) const {
  oops::FieldSet3D fset(validTime_, outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));
  return fset;
}

// -----------------------------------------------------------------------------

void BifourierAnalyticalFilter::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
