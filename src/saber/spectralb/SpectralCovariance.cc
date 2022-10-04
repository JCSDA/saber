/*
 * (C) Crown Copyright 2022 Met Office
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralCovariance.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<SpectralCovariance> makerSpectralCovariance_("spectral covariance");

// -----------------------------------------------------------------------------
SpectralCovariance::SpectralCovariance(const oops::GeometryData & geometryData,
                                       const std::vector<size_t> & variableSizes,
                                       const oops::Variables & inoutVars,
                                       const Parameters_ & params,
                                       const atlas::FieldSet & xb,
                                       const atlas::FieldSet & fg,
                                       const std::vector<atlas::FieldSet> & fsetVec)
  : activeVars_(params.activeVars.value().get_value_or(inoutVars)),
    variance_opt_(params.spectralbParams.value().varianceOpt),
    cs_(), specFunctionSpace_(geometryData.functionSpace())
{
  oops::Log::trace() << classname() << "::SpectralCovariance starting" << std::endl;

  // Initialize CovStat_ErrorCov
  cs_.reset(new CovStat_ErrorCov(variableSizes,
                                 activeVars_,
                                 params.spectralbParams.value()));

  oops::Log::trace() << classname() << "::SpectralCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralCovariance::randomize(atlas::FieldSet & fset) const {
  throw eckit::NotImplemented("SpectralCovariance::randomize: not implemented", Here());
}

// -----------------------------------------------------------------------------

void SpectralCovariance::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  using atlas::array::make_view;
  using atlas::idx_t;

  const idx_t N = specFunctionSpace_.truncation();

  // a spectral convolution with vertical covariances for each total wavenumber
  //   (total wavenumber is n1)
  const atlas::FieldSet & spectralVerticalCovariances = (variance_opt_ ?
                                         cs_->getSpectralVerticalCovariances() :
                                         cs_->getSpectralVerticalCorrelations());
  std::vector<std::string> vertCovNames = spectralVerticalCovariances.field_names();

  const auto zonal_wavenumbers = specFunctionSpace_.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  idx_t i;
  // Only update the fields that were specified in the active variables
  for (const auto & fieldname : activeVars_.variables()) {
    idx_t levels(fieldSet[fieldname].levels());
    auto vertCovView = make_view<const double, 3>(spectralVerticalCovariances[fieldname]);
    auto spfView = make_view<double, 2>(fieldSet[fieldname]);

    i = 0;
    std::vector<double> col(levels), col2(levels);
    double norm;
    for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const idx_t m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
        // note that img stands for imaginary component are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img) {
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[static_cast<std::size_t>(jl)] = spfView(i, jl);
          }
          norm = static_cast<double>((2 * n1 + 1) *
                                     spectralVerticalCovariances[fieldname].shape(0));
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

void SpectralCovariance::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
