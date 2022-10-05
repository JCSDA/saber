/*
 * (C) Crown Copyright 2020-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef SABER_SPECTRALB_SPECTRALBNOINTERP_H_
#define SABER_SPECTRALB_SPECTRALBNOINTERP_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/Partitioner.h"

#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace saber {
namespace spectralb {
namespace detailnointerp {

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid) {
  oops::Log::trace() << "inside createGaussFunctionSpace" << std::endl;
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(1));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
bool createVarianceOpt(const spectralbParameters<MODEL> & params) {
  return (params.varianceOpt.value() != boost::none ?
          params.varianceOpt.value().get() :
          false);
}

}  // namespace detailnointerp
}  // namespace spectralb
}  // namespace saber

// -----------------------------------------------------------------------------
namespace saber {
namespace spectralb {
// -----------------------------------------------------------------------------

template<typename MODEL>
class SpectralBNoInterp {
 public:
  typedef spectralbParameters<MODEL> Parameters_;
  typedef oops::Geometry<MODEL>        Geometry_;
  typedef oops::Increment<MODEL>       Increment_;
  typedef oops::State<MODEL>           State_;

  SpectralBNoInterp(const Geometry_ &,
            const oops::Variables &,
            const Parameters_ &);
  ~SpectralBNoInterp();

  void linearize(const State_ &, const Geometry_ &);
  void multiply(atlas::FieldSet &) const;
  void inverseMultiply(const Increment_ &, Increment_ &) const;
  void randomize(Increment_ &) const;

 private:
  void print(std::ostream &) const;
  std::vector<std::string> vars_;
  std::vector<size_t> varSizes_;
  atlas::StructuredGrid gaussGrid_;
  atlas::functionspace::StructuredColumns gaussFunctionSpace_;
//  atlas::FieldSet  gaussFieldSet_;
  bool variance_opt_;
  std::unique_ptr<const CovStat_ErrorCov<MODEL>> cs_;

  // this method applies the adjoint of the inverse transform
  // then does a convolution with the spectral vertical covariances
  void applySpectralBNoInterp(const atlas::FieldSet &,
                              const atlas::functionspace::Spectral &,
                              const atlas::trans::Trans &,
                              atlas::FieldSet &) const;
};

using atlas::array::make_view;
using atlas::idx_t;

// Note as far as I can see the current implementation of saber blocks assume that the
// geometry object "resol" comes from the model interface.
// However we need to create a Gaussian functionspace here.
// We use the number of vertical levels from "resol"
template<typename MODEL>
SpectralBNoInterp<MODEL>::SpectralBNoInterp(const Geometry_ & resol,
                                            const oops::Variables & vars,
                                            const Parameters_ & params) :
  vars_(vars.variables()),
  varSizes_(resol.variableSizes(vars)),
  gaussGrid_(params.gaussGridUid),
  gaussFunctionSpace_(detailnointerp::createGaussFunctionSpace(gaussGrid_)),
  variance_opt_(detailnointerp::createVarianceOpt(params)),
  cs_(std::make_unique<const CovStat_ErrorCov<MODEL>>(resol, vars, params))
{
  oops::Log::trace() << "SpectralBNoInterp<MODEL>::SpectralBNoInterp done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
SpectralBNoInterp<MODEL>::~SpectralBNoInterp() {
  oops::Log::trace() << "SpectralBNoInterp<MODEL> destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::linearize(const State_ &,
                                 const Geometry_ & resol) {
  oops::Log::trace() << "SpectralBNoInterp<MODEL> linearize" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::multiply(atlas::FieldSet & gaussFieldSet) const {
  oops::Log::trace() << "SpectralBNoInterp<MODEL> multiply start" << std::endl;

  auto N = atlas::GaussianGrid(gaussGrid_).N();

  // assuming that all fields in gaussFieldSet have the same number of levels
  atlas::functionspace::Spectral specFS(
    2*N-1,
    atlas::option::levels(gaussFieldSet[0].levels()));
  atlas::trans::Trans transIFS(gaussFunctionSpace_, specFS);


  // Spectral B
  if (variance_opt_) {
    applySpectralBNoInterp(
      cs_->getSpectralVerticalCovariances(), specFS, transIFS, gaussFieldSet);
  } else {
    applySpectralBNoInterp(
      cs_->getSpectralVerticalCorrelations(), specFS, transIFS, gaussFieldSet);
  }

  gaussFieldSet->haloExchange();

  oops::Log::trace() << "SpectralBNoInterp<MODEL> multiply end"
                     << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::inverseMultiply(const Increment_ & dxin,
                                       Increment_ & dxout) const {
  std::string err_message =
    "saber::SpectralBNoInterp<MODEL>::inverseMultiply not implemented ";
  throw eckit::NotImplemented(err_message, Here());
  oops::Log::trace() << "SpectralBNoInterp<MODEL> inverseMultiply" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::randomize(Increment_ & dx) const {
  oops::Log::trace() << "SpectralBNoInterp<MODEL> randomize" << std::endl;
  std::string err_message =
    "saber::SpectralBNoInterp<MODEL>::randomise not implemented ";
  throw eckit::NotImplemented(err_message, Here());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::print(std::ostream & os) const {
  os << "SpectralBNoInterp<MODEL>::print not implemented";
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SpectralBNoInterp<MODEL>::applySpectralBNoInterp(
    const atlas::FieldSet & spectralVerticalCovariances,
    const atlas::functionspace::Spectral & specFS,
    const atlas::trans::Trans & transIFS,
    atlas::FieldSet & gaussFields) const {
  // the spectral B for each active variable is defined in 3 main steps
  // 1) the adjoint of the inverse spectral transform
  // 2) a spectral convolution with vertical covariances for each total wavenumber
  //   (total wavenumber is n1)
  // 3) the application of the inverse spectral transform

  oops::Log::trace() << "SpectralBNoInterp<MODEL>::applySpectralBNoInterp start" << std::endl;

  std::vector<std::string> fieldNames = gaussFields.field_names();

  idx_t N = specFS.truncation();

  std::vector<std::string> vertCovNames = spectralVerticalCovariances.field_names();

  atlas::FieldSet specFields;

  for (std::size_t f = 0; f < static_cast<std::size_t>(fieldNames.size()); f++) {
    atlas::Field specField =
      specFS.createField<double>(atlas::option::name(fieldNames[f]) |
                                 atlas::option::levels(gaussFields[fieldNames[f]].levels()));
    specFields.add(specField);
  }

  transIFS.invtrans_adj(gaussFields, specFields);

  const auto zonal_wavenumbers = specFS.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  idx_t i;
  for (idx_t f = 0; f < gaussFields.size(); f++) {
    idx_t levels(gaussFields[fieldNames[f]].levels());
    auto vertCovView = make_view<const double, 3>(spectralVerticalCovariances[fieldNames[f]]);
    auto spfView = make_view<double, 2>(specFields[fieldNames[f]]);

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
                                     spectralVerticalCovariances[fieldNames[f]].shape(0));
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

  transIFS.invtrans(specFields, gaussFields);

  oops::Log::trace() << "SpectralBNoInterp<MODEL>::applySpectralBNoInterp end" << std::endl;

  return;
}

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_SPECTRALBNOINTERP_H_
