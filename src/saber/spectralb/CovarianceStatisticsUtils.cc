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

std::vector<std::size_t> getNSpectralBinsFull(const spectralbParameters & params) {
  oops::Variables netCDFVars(params.umatrixNetCDFNames);

  std::vector<std::size_t> nSpectralBinsFull(netCDFVars.size());

  for (std::size_t ivar = 0; ivar < netCDFVars.size(); ++ivar) {
    std::string netCDFVar = netCDFVars[ivar];

    int nbins(0);

    // get the number of spectral bins from the cov file
    covSpectralBins_f90(params.toConfiguration(),
                        static_cast<int>(netCDFVar.size()),
                        netCDFVar.c_str(),
                        nbins);

    nSpectralBinsFull[ivar] = static_cast<std::size_t>(nbins);
  }

  return nSpectralBinsFull;
}

atlas::FieldSet createUMatrices(const oops::Variables & activeVars,
                                const int modelLevels,
                                const std::vector<std::size_t> & nSpectralBinsFull,
                                const spectralbParameters & params) {
  oops::Variables netCDFVars(params.umatrixNetCDFNames);

  atlas::FieldSet spectralUMatrices;

  for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
    std::string var = activeVars[ivar];
    std::string netCDFVar = netCDFVars[ivar];

    auto uMatrix = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nSpectralBinsFull[ivar], modelLevels, modelLevels));

    // vector size
    const int sizeVec = modelLevels * modelLevels * static_cast<int>(nSpectralBinsFull[ivar]);

    std::vector<float> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);

    covSpectralUMatrix_f90(params.toConfiguration(),
                           static_cast<int>(netCDFVar.size()),
                           netCDFVar.c_str(),
                           static_cast<int>(nSpectralBinsFull[ivar]),
                           sizeVec,
                           spectralUMatrix1D[0]);


    auto uMatrixView = atlas::array::make_view<double, 3>(uMatrix);
    std::size_t jn(0);
    for (atlas::idx_t bin = 0; bin < static_cast<atlas::idx_t>(nSpectralBinsFull[ivar]); ++bin) {
      for (atlas::idx_t k1 = 0; k1 < static_cast<atlas::idx_t>(modelLevels); ++k1) {
        for (atlas::idx_t k2 = 0; k2 < static_cast<atlas::idx_t>(modelLevels); ++k2, ++jn) {
          uMatrixView(bin, k1, k2) = spectralUMatrix1D[jn];
        }
      }
    }

    spectralUMatrices.add(uMatrix);
  }

  return spectralUMatrices;
}

// Note - We ideally don't want to have both the Spectral Vertical Covariance and the UMatrices
//        as they double up the memory unnecessarily.
//        There are advantanges to either.
//        Since we are not using a square-root B formulation in JEDI we will not need formally
//        the UMatrices and we can save some unnecessary computation.
//        However, UMatrices are useful for adjoint tests with B
//        where we need the square root of B and for covariance sampling via randomisation.
//        Also it is easier to calculate spectralVerticalCovariances from the UMatrices than
//        vice versa.

atlas::FieldSet createSpectralCovariances(const oops::Variables & activeVars,
                                          const int modelLevels,
                                          const std::vector<std::size_t> & nSpectralBinsFull,
                                          const atlas::FieldSet & spectralUMatrices,
                                          const spectralbParameters & params)
{
  // Read the grid resolution from the grid name.
  std::string gaussGridUid = params.gaussGridUid;
  const auto itFirstDigit = std::find_if(gaussGridUid.begin(), gaussGridUid.end(),
                                         [](char elem){return std::isdigit(elem) != 0;});
  const int nGrid = std::stoi(gaussGridUid.substr(itFirstDigit - gaussGridUid.begin()));
  const int nSpectralBins = 2 * nGrid;

  atlas::FieldSet spectralVerticalCovariances;

  for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
    std::string var = activeVars[ivar];
    ASSERT(static_cast<std::size_t>(nSpectralBins) <= nSpectralBinsFull[ivar]);

    auto spectralVertCov = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nSpectralBins, modelLevels, modelLevels));

    auto spectralVertCovView =
      atlas::array::make_view<double, 3>(spectralVertCov);
    auto uMatrixView = atlas::array::make_view<double, 3>(spectralUMatrices[var]);

    double val;
    for (atlas::idx_t bin = 0; bin < spectralVertCovView.shape(0); ++bin) {
      for (atlas::idx_t k1 = 0; k1 < spectralVertCovView.shape(1); ++k1) {
        for (atlas::idx_t k2 = 0; k2 < spectralVertCovView.shape(2); ++k2) {
          val = 0.0;
          for (atlas::idx_t k3 = 0; k3 < uMatrixView.shape(2); ++k3) {
            val += uMatrixView(bin, k1, k3) * uMatrixView(bin, k2, k3);
          }
          // Crude renormalization assuming the variance is equally distributed across bins.
          spectralVertCovView(bin, k1, k2) = val * nSpectralBins / (nSpectralBinsFull[ivar]);
        }
      }
    }

    spectralVerticalCovariances.add(spectralVertCov);
  }

  return spectralVerticalCovariances;
}

atlas::FieldSet createVerticalSD(const oops::Variables & activeVars,
                                 const int modelLevels,
                                 const atlas::FieldSet & spectralVerticalCovariances) {
  atlas::FieldSet verticalSDs;

  for (std::string var : activeVars.variables()) {
    auto spectralVertCovView =
      atlas::array::make_view<const double, 3>(spectralVerticalCovariances[var]);

    auto verticalSD = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels));

    auto verticalSDView = atlas::array::make_view<double, 1>(verticalSD);

    for (int k = 0; k < modelLevels; ++k) {
      verticalSDView(k) = 0.0;
      for (atlas::idx_t bin = 0; bin < spectralVertCovView.shape(0); ++bin) {
        verticalSDView(k) += spectralVertCovView(bin, k, k);
      }
      verticalSDView(k) = std::sqrt(verticalSDView(k));
    }
    verticalSDs.add(verticalSD);
  }

  return verticalSDs;
}

atlas::FieldSet createSpectralCorrelations(const oops::Variables & activeVars,
                                           const int modelLevels,
                                           const atlas::FieldSet & spectralVerticalCovariances,
                                           const atlas::FieldSet & verticalSDs) {
  atlas::FieldSet spectralCorrelations;

  for (std::string var : activeVars.variables()) {
    auto spectralVertCovView =
      atlas::array::make_view<const double, 3>(spectralVerticalCovariances[var]);

    auto verticalSDView =
      atlas::array::make_view<const double, 1>(verticalSDs[var]);

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

    const atlas::idx_t nSpectralBins = spectralVerticalCovariances[var].shape(0);

    for (int k1 = 0; k1 < modelLevels; ++k1) {
      for (int k2 = 0; k2 < modelLevels; ++k2) {
        correlationScalingView(k1, k2) =
          static_cast<double>(nSpectralBins) /
          (verticalSDView(k1) * verticalSDView(k2));
      }
    }

    for (atlas::idx_t bin = 0; bin < spectralVertCorrel.shape(0); ++bin) {
      for (atlas::idx_t k1 = 0; k1 < spectralVertCorrel.shape(1); ++k1) {
        for (atlas::idx_t k2 = 0; k2 < spectralVertCorrel.shape(2); ++k2) {
          spectralVertCorrelView(bin, k1, k2) = spectralVertCovView(bin, k1, k2) *
            correlationScalingView(k1, k2);
        }
      }
    }

    spectralCorrelations.add(spectralVertCorrel);
  }

  return spectralCorrelations;
}


}  // namespace spectralb
}  // namespace saber
