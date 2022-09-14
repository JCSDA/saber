/*
 * (C) Crown Copyright 2020-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef SABER_SPECTRALB_SPECTRALB_H_
#define SABER_SPECTRALB_SPECTRALB_H_

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

#include "atlas/redistribution/Redistribution.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/Trans.h"

#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/spectralb/CovarianceStatistics.h"
#include "saber/spectralb/spectralbParameters.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace saber {
namespace spectralb {
namespace detail {

atlas::functionspace::StructuredColumns
    createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid) {
  oops::Log::trace() << "inside createGaussFunctionSpace" << std::endl;
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(1));
}

std::shared_ptr<atlas::FieldSet> allocateGaussFieldset(
    const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
    const std::vector<std::string> & gaussNames,
    const std::shared_ptr<const atlas::FieldSet> & fieldsetIn) {

  oops::Log::trace() << "allocateGaussFieldset starting" << std::endl;

  // create gauss FieldSet (fields are zeroed)
  auto gaussFieldSet = std::make_shared<atlas::FieldSet>();
  std::vector<std::string> varin = (*fieldsetIn).field_names();

  for (std::size_t i = 0; i < varin.size(); ++i) {
    atlas::Field gaussField =
      gaussFunctionSpace.createField<double>(
        atlas::option::name(gaussNames[i]) |
        atlas::option::levels((*fieldsetIn)[varin[i]].levels()) );

    gaussField.haloExchange();

    (*gaussFieldSet).add(gaussField);
  }

  gaussFieldSet->haloExchange();
  oops::Log::trace() << "allocateGaussFieldset done" << std::endl;
  return gaussFieldSet;
}

// -----------------------------------------------------------------------------

atlas::Grid createOutputGrid(const spectralbParameters & params) {
  std::string gridName((params.outputGridUid.value() != boost::none ?
                        params.outputGridUid.value().get() :
                        params.gaussGridUid));
  return atlas::Grid(gridName);
}

// -----------------------------------------------------------------------------

atlas::FunctionSpace createOutputFunctionSpace(const atlas::FieldSet & fset) {
  return fset[0].functionspace();
}

// -----------------------------------------------------------------------------

bool createVarianceOpt(const spectralbParameters & params) {
  return (params.varianceOpt.value() != boost::none ?
          params.varianceOpt.value().get() :
          false);
}

}  // namespace detail
}  // namespace spectralb
}  // namespace saber

// -----------------------------------------------------------------------------
namespace saber {
namespace spectralb {
// -----------------------------------------------------------------------------

class SpectralB {
 public:
  typedef spectralbParameters Parameters_;

  SpectralB(const atlas::FunctionSpace &,
            const std::vector<size_t> &,
            const oops::Variables &,
            const Parameters_ &);
  ~SpectralB();

  void multiply_InterpAndCov(atlas::FieldSet &) const;
  static atlas::FieldSet createFieldsSpace(const atlas::FunctionSpace &,
                                           const std::vector<size_t> &,
                                           const oops::Variables &);

 private:
  std::shared_ptr<const atlas::FieldSet> modelFieldSet_;
  std::vector<std::string> gaussNames_;
  atlas::StructuredGrid gaussGrid_;
  atlas::functionspace::StructuredColumns gaussFunctionSpace_;
  std::shared_ptr<atlas::FieldSet> gaussFieldSet_;
  saber::interpolation::AtlasInterpWrapper interp_;
  bool variance_opt_;
  std::unique_ptr<const CovStat_ErrorCov> cs_;

  // this method applies the adjoint of the inverse transform
  // then does a convolution with the spectral vertical covariances
  void applySpectralB(const atlas::FieldSet &,
                      const atlas::functionspace::Spectral &,
                      const atlas::trans::Trans &,
                      atlas::FieldSet &) const;
};

using atlas::array::make_view;
using atlas::idx_t;

SpectralB::SpectralB(const atlas::FunctionSpace & functionSpace,
                     const std::vector<size_t> & variableSizes,
                     const oops::Variables & vars,
                     const Parameters_ & params) :
  modelFieldSet_(std::make_shared<const atlas::FieldSet>(createFieldsSpace(functionSpace,
                                                                           variableSizes,
                                                                           vars))),
  gaussNames_(vars.variables()),
  gaussGrid_(params.gaussGridUid),
  gaussFunctionSpace_(detail::createGaussFunctionSpace(gaussGrid_)),
  gaussFieldSet_(detail::allocateGaussFieldset(gaussFunctionSpace_, gaussNames_, modelFieldSet_)),
  interp_(atlas::grid::Partitioner(new TransPartitioner()), gaussFunctionSpace_,
    detail::createOutputGrid(params), detail::createOutputFunctionSpace(*modelFieldSet_)),
  variance_opt_(detail::createVarianceOpt(params)),
  cs_(std::make_unique<const CovStat_ErrorCov>(variableSizes, vars, params))
{
  oops::Log::trace() << "SpectralB::SpectralB done" << std::endl;
}

// -----------------------------------------------------------------------------

SpectralB::~SpectralB() {
  oops::Log::trace() << "SpectralB destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralB::multiply_InterpAndCov(atlas::FieldSet & modelGridFieldSet) const {
  oops::Log::trace() << "SpectralB multiply_InterpAndCov start" << std::endl;

  auto N = atlas::GaussianGrid(gaussGrid_).N();

  // assuming that all fields in modelGridFieldSet have the same number of levels
  atlas::functionspace::Spectral specFS(2*N-1,
                                        atlas::option::levels(modelGridFieldSet[0].levels()));
  atlas::trans::Trans transIFS(gaussFunctionSpace_, specFS);

  interp_.executeAdjoint(*gaussFieldSet_, modelGridFieldSet);

  // Spectral B
  if (variance_opt_) {
    applySpectralB(cs_->getSpectralVerticalCovariances(), specFS, transIFS, *gaussFieldSet_);
  } else {
    applySpectralB(cs_->getSpectralVerticalCorrelations(), specFS, transIFS, *gaussFieldSet_);
  }

  gaussFieldSet_->haloExchange();

  interp_.execute(*gaussFieldSet_, modelGridFieldSet);

  oops::Log::trace() << "SpectralB multiply_InterpAndCov end"
                     << std::endl;
}

// -----------------------------------------------------------------------------

atlas::FieldSet SpectralB::createFieldsSpace(const atlas::FunctionSpace & functionSpace,
                                             const std::vector<size_t> & variableSizes,
                                             const oops::Variables & vars) {
  oops::Log::trace() << "start createFieldsSpace:: Rank, Variables = "
                     << atlas::mpi::rank() << " " << vars << std::endl;

  atlas::FieldSet fields;

  ASSERT(vars.size() > 0);

  std::vector<size_t> sizes = variableSizes;

  for (unsigned int i = 0; i < vars.size(); ++i) {
    atlas::FunctionSpace nodesFs = functionSpace;

    atlas::Field tempField =
      nodesFs.createField<double>(atlas::option::name(vars[i]) |
                                  atlas::option::levels(sizes[i]));

    auto view = atlas::array::make_view<double, 2>(tempField);
    view.assign(0.0);

    fields.add(tempField);
  }
  return fields;
}

// -----------------------------------------------------------------------------

void SpectralB::applySpectralB(
    const atlas::FieldSet & spectralVerticalCovariances,
    const atlas::functionspace::Spectral & specFS,
    const atlas::trans::Trans & transIFS,
    atlas::FieldSet & gaussFields) const {
  // the spectral B for each active variable is defined in 3 main steps
  // 1) the adjoint of the inverse spectral transform
  // 2) a spectral convolution with vertical covariances for each total wavenumber
  //   (total wavenumber is n1)
  // 3) the application of the inverse spectral transform

  oops::Log::trace() << "SpectralB::applySpectralB start" << std::endl;

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
  const int nb_zonal_wavenumbers = zonal_wavenumbers.size();

  int i;
  for (idx_t f = 0; f < gaussFields.size(); f++) {
    int levels(gaussFields[fieldNames[f]].levels());
    auto vertCovView = make_view<const double, 3>(spectralVerticalCovariances[fieldNames[f]]);
    auto spfView = make_view<double, 2>(specFields[fieldNames[f]]);

    i = 0;
    std::vector<double> col(levels), col2(levels);
    double norm;
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(N); ++n1) {
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

  oops::Log::trace() << "SpectralB::applySpectralB end" << std::endl;

  return;
}

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_SPECTRALB_H_
