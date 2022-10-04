/*
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "saber/spectralb/CovarianceStatisticsUtils.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "saber/spectralb/spectralb_covstats_interface.h"
#include "saber/spectralb/spectralbParameters.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace saber {
namespace spectralb {

std::vector<std::size_t> getNetCDFSpectralBins(const spectralbParameters & params) {
  oops::Variables netCDFVars(params.umatrixNetCDFNames);

  std::vector<std::size_t> netCDFSpectralBins(netCDFVars.size());

  for (std::size_t i = 0; i < netCDFVars.size(); ++i) {
    std::string netCDFVar = netCDFVars[i];

    int nbins(0);

    // get the number of spectral bins from the cov file
    covSpectralBins_f90(params.toConfiguration(),
                        static_cast<int>(netCDFVar.size()),
                        netCDFVar.c_str(),
                        nbins);

    netCDFSpectralBins[i] = static_cast<std::size_t>(nbins);
  }

  return netCDFSpectralBins;
}

atlas::FieldSet createUMatrices(const oops::Variables & activeVars,
                                const int modelLevels,
                                const std::vector<std::size_t> & netCDFSpectralBins,
                                const spectralbParameters & params) {
  // we are currently not using globalNons_ here
  // however in the near future we may want to scale the spectral covariances
  // so that they can be run at lower
  // spectral resolutions.  That is the reason for the globalNLons argument.
  oops::Variables netCDFVars(params.umatrixNetCDFNames);

  atlas::FieldSet spectralUMatrices;

  for (std::size_t i = 0; i < activeVars.size(); ++i) {
    std::string var = activeVars[i];
    std::string netCDFVar = netCDFVars[i];

    auto field = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(netCDFSpectralBins[i], modelLevels, modelLevels));

    // vector size
    int sizeVec = modelLevels * modelLevels * static_cast<int>(netCDFSpectralBins[i]);

    std::vector<float> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);

    covSpectralUMatrix_f90(params.toConfiguration(),
                           static_cast<int>(netCDFVar.size()),
                           netCDFVar.c_str(),
                           static_cast<int>(netCDFSpectralBins[i]),
                           sizeVec,
                           spectralUMatrix1D[0]);


    auto fview = atlas::array::make_view<double, 3>(field);
    std::size_t jn(0);
    for (atlas::idx_t b = 0; b < static_cast<atlas::idx_t>(netCDFSpectralBins[i]); ++b) {
      for (atlas::idx_t j = 0; j < static_cast<atlas::idx_t>(modelLevels); ++j) {
        for (atlas::idx_t k = 0; k < static_cast<atlas::idx_t>(modelLevels); ++k, ++jn) {
          fview(b, j, k) = spectralUMatrix1D[jn];
        }
      }
    }

    spectralUMatrices.add(field);
  }

  return spectralUMatrices;
}

// Note - we ideally don't want to have both the Spectral Vertical Covariance and the UMatrices
//        as they double up the memory unnecessarily.
//        There are advantanges to either.
//        Since we are not using a square-root B formulation in JEDI we will not need formally
//        the UMatrices and we can save some unnecessary computation.
//        However, UMatrices are useful for adjoint tests with B
//        where we need the square root of B and for covariance sampling via randomsiation.
//        Also it is easier the calculate spectralVerticalCovariances from the UMatrices than
//        visa versa.
atlas::FieldSet createSpectralCovariances(const oops::Variables & activeVars,
                                          const int modelLevels,
                                          const std::vector<std::size_t> & netCDFSpectralBins,
                                          const atlas::FieldSet & spectralUMatrices,
                                          const spectralbParameters & params)
{
  std::string gaussGridUid = params.gaussGridUid;
  // extract resolution find position where string has a number
  std::size_t j(0);
  for (std::size_t i = 0; i < gaussGridUid.size(); ++i) {
    if (std::isdigit(gaussGridUid[i]) == 0) {
      ++j;
    }
  }
  int nBins = 2 * std::stoi(gaussGridUid.substr(j, gaussGridUid.size() - j));

  atlas::FieldSet spectralVerticalCovariances;

  for (std::size_t i = 0; i < activeVars.size(); ++i) {
    std::string var = activeVars[i];

    auto spectralVertCov = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nBins, modelLevels, modelLevels));

    auto spectralVertCovView =
      atlas::array::make_view<double, 3>(spectralVertCov);
    auto uMatrixView = atlas::array::make_view<double, 3>(spectralUMatrices[var]);

    double val;
    for (atlas::idx_t b = 0; b < spectralVertCovView.shape(0); ++b) {
      for (atlas::idx_t r = 0; r < spectralVertCovView.shape(1); ++r) {
        for (atlas::idx_t c = 0; c < spectralVertCovView.shape(2); ++c) {
          val = 0.0;
          for (atlas::idx_t s = 0; s < spectralVertCovView.shape(2); ++s) {
            val += uMatrixView(b, r, s) * uMatrixView(b, c, s);
          }
          spectralVertCovView(b, r, c) = val * nBins / (netCDFSpectralBins[i]);
        }
      }
    }

    spectralVerticalCovariances.add(spectralVertCov);
  }

  return spectralVerticalCovariances;
}

atlas::FieldSet createSpectralSD(const oops::Variables & activeVars,
                                 const int modelLevels,
                                 const atlas::FieldSet & spectralVerticalCovariances) {
  atlas::FieldSet spectralSDs;

  for (std::string var : activeVars.variables()) {
    auto spectralVerticalCovarianceView =
      atlas::array::make_view<const double, 3>(spectralVerticalCovariances[var]);

    auto spectralSD = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels));

    auto spectralSDView = atlas::array::make_view<double, 1>(spectralSD);

    for (int k = 0; k < modelLevels; ++k) {
      spectralSDView(k) = 0.0;
      for (atlas::idx_t b = 0; b < spectralVerticalCovarianceView.shape(0); ++b) {
        spectralSDView(k) += spectralVerticalCovarianceView(b, k, k);
      }
      spectralSDView(k) = std::sqrt(spectralSDView(k));
    }
    spectralSDs.add(spectralSD);
  }

  return spectralSDs;
}

atlas::FieldSet createSpectralCorrelations(const oops::Variables & activeVars,
                                           const int modelLevels,
                                           const atlas::FieldSet & spectralVerticalCovariances,
                                           const atlas::FieldSet & spectralSDs) {
  atlas::FieldSet spectralCorrelations;

  for (std::string var : activeVars.variables()) {
    auto spectralVertCovView =
      atlas::array::make_view<const double, 3>(spectralVerticalCovariances[var]);

    auto spectralSDView =
      atlas::array::make_view<const double, 1>(spectralSDs[var]);

    auto spectralVertCorrel = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(spectralVertCovView.shape(0),
                               spectralVertCovView.shape(1),
                               spectralVertCovView.shape(2)));

    auto correlationScaling = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels, modelLevels));

    auto spectralVertCorrelView =
      atlas::array::make_view<double, 3>(spectralVertCorrel);

    auto correlationScalingView =
      atlas::array::make_view<double, 2>(correlationScaling);

    atlas::idx_t nBins = spectralVerticalCovariances[var].shape(0);

    for (int k1 = 0; k1 < modelLevels; ++k1) {
      for (int k2 = 0; k2 < modelLevels; ++k2) {
        correlationScalingView(k1, k2) =
          static_cast<double>(nBins)/
          (spectralSDView(k1) * spectralSDView(k2));
      }
    }

    for (atlas::idx_t b = 0; b < spectralVertCorrel.shape(0); ++b) {
      for (atlas::idx_t r = 0; r < spectralVertCorrel.shape(1); ++r) {
        for (atlas::idx_t c = 0; c < spectralVertCorrel.shape(2); ++c) {
          spectralVertCorrelView(b, r, c) = spectralVertCovView(b, r, c) *
            correlationScalingView(r, c);
        }
      }
    }

    spectralCorrelations.add(spectralVertCorrel);
  }

  return spectralCorrelations;
}


}  // namespace spectralb
}  // namespace saber
