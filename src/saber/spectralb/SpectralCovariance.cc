/*
 * (C) Crown Copyright 2022 Met Office
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

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<SpectralCovariance> makerSpectralCovariance_("spectral covariance");

// -----------------------------------------------------------------------------
SpectralCovariance::SpectralCovariance(const oops::GeometryData & geometryData,
                                       const std::vector<size_t> & variableSizes,
                                       const oops::Variables & centralVars,
                                       const eckit::Configuration & covarConf,
                                       const Parameters_ & params,
                                       const atlas::FieldSet & xb,
                                       const atlas::FieldSet & fg,
                                       const util::DateTime & validTimeOfXbFg,
                                       const size_t & timeRank)
  : params_(params),
    variableSizes_(variableSizes),
    activeVars_(params.activeVars.value().get_value_or(centralVars)),
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

void SpectralCovariance::testUUtConsistency(const double & consistencyTolerance,
                                            const double & adjointTolerance) const {
  // Test that UU^t = B and that U and U^t are adjoint.
  oops::Log::trace() << "SpectralCovariance::testUUtConsistency starting" << std::endl;

  // Step 1/2: Test that UU^t = B
  // Create random FieldSet y
  atlas::FieldSet y =  util::createRandomFieldSet(geometryData_.comm(),
                                                  geometryData_.functionSpace(),
                                                  variableSizes_,
                                                  activeVars_.variables(),
                                                  timeRank_);

  // Apply forward multiplication using B or UU^t
  atlas::FieldSet By = util::copyFieldSet(y);
  this->multiply(By);

  atlas::FieldSet Uty = util::copyFieldSet(y);
  this->multiplyUMatrixAD(Uty);

  atlas::FieldSet UUty = util::copyFieldSet(Uty);
  this->multiplyUMatrix(UUty);

  // Check that B y == UU^t y
  oops::Log::test() << "UU^t consistency test";
  const oops::Variables variables1(By.field_names());
  const oops::Variables variables2(UUty.field_names());
  if (variables1 != variables2) {
      oops::Log::test() << " failed" << std::endl;
      ABORT("UU^t consistency test failure, wrong variables");
  }

  for (const atlas::Field & field1 : By) {
    // Get associated field
    const atlas::Field & field2 = UUty.field(field1.name());
    // Compare ranks
    if (field1.rank() != field2.rank()) {
        oops::Log::test() << " failed" << std::endl;
        ABORT("UU^t consistency test failure, wrong ranks");
    }
    // Compare data
    if (field1.rank() == 2) {
      const auto view1 = atlas::array::make_view<double, 2>(field1);
      const auto view2 = atlas::array::make_view<double, 2>(field2);
      for (int jnode = 0; jnode < field1.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field1.shape(1); ++jlevel) {
          const double d1 = view1(jnode, jlevel);
          const double d2 = view2(jnode, jlevel);
          if (0.5*abs(d1-d2)/(d1+d2) >= consistencyTolerance) {
            oops::Log::test() << " failed" << std::endl;
            ABORT("UU^t consistency test failure, differing data");
          }
        }
      }
    } else {
      oops::Log::test() << " failed" << std::endl;
      ABORT("UU^t consistency test failure, expected rank 2");
    }
  }

  oops::Log::test() << " passed" << std::endl;

  // Step 2/2: Adjoint test for U and U^t
  // Generate random field set x in spectral space
  atlas::FieldSet x = util::createRandomFieldSet(geometryData_.comm(),
                                                 geometryData_.functionSpace(),
                                                 variableSizes_,
                                                 activeVars_.variables());
  // Apply U
  atlas::FieldSet Ux = util::copyFieldSet(x);
  this->multiplyUMatrix(Ux);

  // Compute adjoint test
  const double dp1 = util::dotProductFieldSets(y, Ux, activeVars_.variables(),
                                               geometryData_.comm());
  const double dp2 = util::dotProductFieldSets(Uty, x, activeVars_.variables(),
                                               geometryData_.comm());
  oops::Log::info() << "Info     : Adjoint test: <y ; Ux> = " << dp1
                    << ": <U^t y ; x> = " << dp2 << std::endl;
  oops::Log::test() << "Adjoint test for U and U^t";
  if (0.5*abs(dp1-dp2)/(dp1+dp2) < adjointTolerance) {
    oops::Log::test() << " passed" << std::endl;
  } else {
    oops::Log::test() << " failed" << std::endl;
    ABORT("Adjoint test failure");
  }

  oops::Log::trace() << "SpectralCovariance::testUUtConsistency done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::randomize(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Overwrite input fieldSet with random numbers in spectral space
  atlas::FieldSet newFieldSet = util::createRandomFieldSet(geometryData_.comm(),
                                                           geometryData_.functionSpace(),
                                                           variableSizes_,
                                                           activeVars_.variables());
  for (auto & var : activeVars_.variables()) {
    fieldSet[var] = newFieldSet[var];
  }

  // Apply U matrix
  multiplyUMatrix(fieldSet);

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

  idx_t i;
  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t levels(fieldSet[var].levels());
    auto vertCovView = make_view<const double, 3>(spectralVerticalCovariances[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);

    i = 0;
    std::vector<double> col(levels), col2(levels);
    double norm;
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
          norm = static_cast<double>((2 * n1 + 1) * vertCovView.shape(0));
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

void SpectralCovariance::multiplyUMatrix(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyUMatrix starting" << std::endl;
  using atlas::array::make_view;
  using atlas::idx_t;

  const idx_t N = specFunctionSpace_.truncation();

  // 1D convolution with spectral vertical covariances for each total wavenumber n1
  const atlas::FieldSet & spectralUMatrices = (variance_opt_ ?
                                         cs_->getSpectralUMatrices() :
                                         cs_->getSpectralCorrelUMatrices());

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  idx_t i;
  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t levels(fieldSet[var].levels());
    auto UMatrixView = make_view<const double, 3>(spectralUMatrices[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);
    const int nSpectralBinsFull = spectralUMatrices[var].shape(0);

    i = 0;
    std::vector<double> col(levels), col2(levels);
    double norm;
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
          norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[static_cast<std::size_t>(jl)] / norm;
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplyUMatrix done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::multiplyUMatrixAD(atlas::FieldSet & fieldSet) const {
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

  idx_t i;
  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars_.variables()) {
    idx_t levels(fieldSet[var].levels());
    auto UMatrixView = make_view<const double, 3>(spectralUMatrices[var]);
    auto spfView = make_view<double, 2>(fieldSet[var]);
    const int nSpectralBinsFull = spectralUMatrices[var].shape(0);

    i = 0;
    std::vector<double> col(levels), col2(levels);
    double norm;
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
          norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
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

void SpectralCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Initialize CovStat_ErrorCov
  cs_.reset(new CovStat_ErrorCov(variableSizes_,
                                 activeVars_,
                                 *params_.readParams.value()));

  // Check consistency of UU^t and B
  static bool uut_first_pass = true;
  if (params_.readParams.value()->uutConsistencyTest && uut_first_pass) {
      testUUtConsistency(params_.readParams.value()->consistencyTolerance);
      uut_first_pass = false;
  }
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
