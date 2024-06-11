/*
 * (C) Crown Copyright 2022-2024 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "saber/vader/CovarianceStatisticsUtils.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/CoordinateEnums.h"

#include "mo/constants.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "saber/spectralb/spectralb_covstats_interface.h"
#include "saber/vader/movader_covstats_interface.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

/// \details This extracts the hp_gp regression matrix for a number of
///          of overlapping latitude bands from the operational covariance
///          statistics file. The matrices are stored as a single field.
///
/// B = (vertical regression matrix bin_0)
///     (vertical regression matrix bin_1)
///     (          ...                   )
///     (vertical regression matrix bin_m)
/// Since each matrix is square we can easily infer the bin index from the row index
/// First index of vertRegView is bin_index * number of levels + level index,
///     the second is number of levels associated with matrix column.
///
/// The interpolation weights are calculated for each grid point location
/// The second index relates to the bin index.
/// We ensure that across all bins for a grid point we sum to 1.

atlas::Field createGpRegressionMatrices(const std::string & covFileName,
                                        const std::size_t gpBins,
                                        const std::size_t modelLevels) {
  auto regMatrixFld = atlas::Field(std::string("vertical_regression_matrices"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(gpBins * modelLevels, modelLevels));

  std::vector<float> regresssionMatrices1D(modelLevels * modelLevels * gpBins, 0.0);

  covRegressionMatrices_f90(covFileName.size(),
                            covFileName.c_str(),
                            static_cast<int>(modelLevels),
                            static_cast<int>(gpBins),
                            static_cast<int>(regresssionMatrices1D.size()),
                            regresssionMatrices1D[0]);

  auto fview = atlas::array::make_view<double, 2>(regMatrixFld);
  std::size_t jn(0);
  for (atlas::idx_t j = 0; j < static_cast<atlas::idx_t>(modelLevels * gpBins); ++j) {
    for (atlas::idx_t k = 0; k < static_cast<atlas::idx_t>(modelLevels); ++k, ++jn) {
      fview(j, k) = regresssionMatrices1D[jn];
    }
  }

  return regMatrixFld;
}

// -----------------------------------------------------------------------------

std::vector<double> interpWeights(std::vector<std::vector<double>> & regWeights,
                                  std::vector<std::vector<double>> & latValues,
                                  double latDest) {
  std::vector<double> interpWeight(regWeights.size(), 0.0);
  float invCovLatDelta =
    1.0/(latValues[0][1] - latValues[0][0]);

  // note that there are some similarites between this code and
  // what is in SVP.cc. At some point it might make sense to convert this
  // into a generic interpolator
  auto inRange = [&](const double latPt, std::size_t b, std::size_t bRows) {
    return latPt >= latValues[b][0] &&
           latPt <= latValues[b][bRows-1];
  };

  auto normalisedLat = [&](const double latPt, std::size_t b) {
    return (latPt - latValues[b][0]) * invCovLatDelta;
  };

  auto index = [](const double normalisedLat, std::size_t bRows) {
    std::size_t i = static_cast<std::size_t>(normalisedLat);
    return (i == bRows - 1 ? bRows - 2 : i);
  };

  auto weight = [index](const double normalisedLat, std::size_t bRows) {
      std::size_t i = index(normalisedLat, bRows);
      return (normalisedLat - static_cast<float>(i));
  };

  // interp lambda function
  auto interp = [](const double weight, const double T1, const double T2) {
    return (weight * T2 + (1 - weight) * T1);
  };

  for (std::size_t b = 0; b < regWeights.size(); b++) {
    std::size_t bRows = static_cast<std::size_t>(regWeights[b].size());
    if (inRange(latDest, b, bRows)) {
      float normLat = normalisedLat(latDest, b);
      std::size_t indx = index(normLat, bRows);
      float w = weight(normLat, bRows);
      interpWeight[b] = interp(w, regWeights[b][indx], regWeights[b][indx+1]);
    }
  }

  return interpWeight;
}

// -----------------------------------------------------------------------------

atlas::Field createGpRegressionWeights(const atlas::FunctionSpace & functionSpace,
                                       const atlas::FieldSet & fields,
                                       const std::string & covFileName,
                                       const std::size_t covGlobalNLats,
                                       const std::size_t gpBins) {
  // get interpolation weights from file
  atlas::idx_t horizPts = fields[0].shape(0);

  std::vector<float> regresssionWeights1D(covGlobalNLats * gpBins, 0.0);
  std::vector<float> covLatitudesVec(covGlobalNLats);
  std::vector<int> startVec(gpBins, 0);
  std::vector<int> lenVec(gpBins, 0);

  covRegressionWeights_f90(covFileName.size(),
                           covFileName.c_str(),
                           static_cast<int>(covGlobalNLats),
                           static_cast<int>(gpBins),
                           static_cast<int>(regresssionWeights1D.size()),
                           startVec[0],
                           lenVec[0],
                           covLatitudesVec[0],
                           regresssionWeights1D[0]);

  std::vector<std::vector<double>> latValues;
  std::vector<std::vector<double>> regWeights;

  std::size_t tot(0);
  for (std::size_t b = 0; b < gpBins; ++b) {
  // find number of non-zero weights for each bin.

    std::vector<double> latValuesBin(lenVec[b]);
    std::vector<double> regWeightsBin(lenVec[b]);

    for (int i = 0; i < lenVec[b]; ++i) {
      latValuesBin[i] = covLatitudesVec[startVec[b] + i];
      regWeightsBin[i] = regresssionWeights1D[tot + i];
    }

    latValues.push_back(latValuesBin);
    regWeights.push_back(regWeightsBin);

    tot += static_cast<std::size_t>(lenVec[b]);
  }

  // need to look over horiz latitude points to calculate gp regression.
  // the horizontal points need to be PE decomposed.
  auto interWgtFld =
    functionSpace.createField<double>(atlas::option::name("interpolation_weights") |
                                      atlas::option::levels(gpBins));

  auto lonlatView = atlas::array::make_view<double, 2>(functionSpace.lonlat());
  auto interWgtFldView = atlas::array::make_view<double, 2>(interWgtFld);

  std::vector<double> tempWgt(gpBins);
  for (atlas::idx_t h = 0; h < horizPts; ++h) {
    tempWgt = interpWeights(regWeights, latValues, lonlatView(h, atlas::LAT));

    double invWeightTot = 1.0 /
          (std::accumulate(begin(tempWgt), end(tempWgt), 0.0));

    for (std::size_t b = 0; b < gpBins; ++b) {
      interWgtFldView(h, b) = tempWgt[b] * invWeightTot;
    }
  }

  return interWgtFld;
}

// -----------------------------------------------------------------------------

void interpMuStats(atlas::FieldSet & augmentedStateFieldSet,
                   const atlas::Field & covFld) {
  // variable name of field to be populated
  // (created by removing the "stats" from the end of string)
  std::string covFldName(covFld.name());
  std::string varName = covFldName.substr(0, covFldName.size() - 5);

  atlas::idx_t muBins = covFld.shape()[1];

  double inverseBinDelta = 1.0 / mo::constants::rHTBin;

  // note that there are some similarites between this code and
  // what is in SVP.cc. At some point it might make sense to convert this
  // into a generic interpolator

  auto normalisedField = [&](const double rhVal) {
    return (rhVal - mo::constants::MinRhRef ) * inverseBinDelta;
  };

  auto index = [](const double normalisedField, const int muBins) {
    int i = static_cast<int>(normalisedField);
    // Restrict i to lie within [0, muBins-2] inclusive
    return std::max(std::min(i, muBins - 2), 0);
  };

  auto weight = [](const double normalisedField, const int index) {
    return (normalisedField - static_cast<double>(index));
  };

  // interp lambda function
  auto interp = [](const double weight, const double beforePt, const double afterPt) {
    return ( weight * afterPt + (1 - weight) * beforePt);
  };

  auto RHtView = atlas::array::make_view<double, 2>(augmentedStateFieldSet["rht"]);
  auto muView = atlas::array::make_view<double, 2>(augmentedStateFieldSet[varName]);
  auto covFldView = atlas::array::make_view<double, 2>(covFld);

  for (atlas::idx_t jn = 0; jn < augmentedStateFieldSet["rht"].shape(0); ++jn) {
    for (int jl = 0; jl < augmentedStateFieldSet["rht"].shape(1); ++jl) {
      double normField = normalisedField(static_cast<double>(RHtView(jn, jl)));
      int indx = index(normField, muBins);
      double w = weight(normField, indx);
      muView(jn, jl) = static_cast<double>(interp(w, covFldView(jl, indx),
                                                     covFldView(jl, indx + 1)));
    }
  }
  if (varName.compare("muA") == 0) {
    for (atlas::idx_t jn = 0; jn < augmentedStateFieldSet[varName].shape(0); ++jn) {
      for (int jl = 0; jl < augmentedStateFieldSet[varName].shape(1); ++jl) {
        // note -  the actual mu inverse field in VAR involves an additional step
        //         that involves interpolation of a "mu_table".
        //         What is here is a simplification.
        muView(jn, jl) = 1.0 / muView(jn, jl);
      }
    }
  }
}


void populateMuA(atlas::FieldSet & augmentedStateFieldSet,
                 const atlas::Field & covFld) {
  atlas::idx_t modelLevels(covFld.shape(0));
  atlas::idx_t muBins(covFld.shape(1));
  auto covFldView = atlas::array::make_view<double, 2>(covFld);

  auto invStatsFld = atlas::Field("invMuStats",
                                  atlas::array::make_datatype<double>(),
    atlas::array::make_shape(modelLevels, muBins + 1));
  auto invStatsFldView = atlas::array::make_view<double, 2>(invStatsFld);

  std::vector<double> logStats(muBins, 0.0);
  std::vector<double> sigmaExtended(muBins + 2, 0.0);
  for (atlas::idx_t jl = 0; jl < modelLevels; ++jl) {
    // calculate the natural logarithm of the original normalisation statistics.
    for (atlas::idx_t bin = 0; bin < muBins; ++bin) {
      logStats[bin] = std::log(covFldView(jl, bin));
    }
    // The bottom bin is often poorly sampled, so extrapolate down from bin 2
    sigmaExtended[0] = covFldView(jl, 1) * mo::constants::effectiveRNegative /
       (mo::constants::MinRhRef + 1.5 * mo::constants::rHTBin);

    for (std::size_t ja = 1; ja < sigmaExtended.size(); ++ja) {
      const double rhTColumn =  (static_cast<double>(ja)-0.5) * mo::constants::rHTBin;
      if (rhTColumn <= mo::constants::MinRhRef) {
        // make sigma tend to zero proportionally to rhTColumn
        // The bottom bin is often poorly sampled, so extrapolate down from bin 2
        sigmaExtended[ja] = covFldView(jl, 1) * rhTColumn /
          (mo::constants::MinRhRef + 1.5 * mo::constants::rHTBin);
      } else {
        // interpolate or extrapolate log(sigma) in order to ensure sigma>0
        double s = (rhTColumn - mo::constants::MinRhRef) / mo::constants::rHTBin + 0.5;
        const int j = std::max(1, std::min(muBins-1, static_cast<int>(s)));
        s -= j;
        sigmaExtended[ja] = std::exp(logStats[j-1] + s * (logStats[j] - logStats[j-1]));
      }
    }

    const double tol = mo::constants::TolMonotonicity * mo::constants::rHTBin;
    for (int jb = 0; jb < muBins + 1; ++jb) {
      const double s = (sigmaExtended[jb] + sigmaExtended[jb+1])/2.0;
      invStatsFldView(jl, jb) = - mo::constants::rHTBin / s;
      invStatsFldView(jl, jb)  =
        invStatsFldView(jl, jb)  > - tol ? -tol : invStatsFldView(jl, jb);
    }
  }

  // interpolate invStats onto model mesh
  auto RHtView = atlas::array::make_view<double, 2>(augmentedStateFieldSet["rht"]);
  auto muView = atlas::array::make_view<double, 2>(augmentedStateFieldSet["muA"]);

  const double rhTlow = -0.5 * mo::constants::rHTBin;
  const double rhThigh = rhTlow + (muBins + 1) *  mo::constants::rHTBin;

  for (atlas::idx_t jl = 0; jl < modelLevels; ++jl) {
    for (atlas::idx_t jn = 0; jn < RHtView.shape(0); ++jn) {
      double wgt = (std::max(rhTlow, std::min(rhThigh, RHtView(jn, jl))) - rhTlow) /
                    mo::constants::rHTBin;
      int indx = std::min(static_cast<int>(wgt), static_cast<int>(invStatsFldView.shape(1)) - 1);
      wgt -= static_cast<double>(indx);

      if (indx == 0) {
        muView(jn, jl) = - invStatsFldView(jl, 0) / mo::constants::rHTBin;
      } else {
        muView(jn, jl) = ((1.0 - wgt) *
                          (- invStatsFldView(jl, indx-1)) +
                          wgt * (- invStatsFldView(jl, indx))) / mo::constants::rHTBin;
      }
    }
  }
}

// -----------------------------------------------------------------------------

/// \details This extracts the hp_gp regression matrix for a number of
///          of overlapping latitude bands from the operational covariance
///          statistics file. The matrices are stored as a single field.
///
/// B = (vertical regression matrix bin_0)
///     (vertical regression matrix bin_1)
///     (          ...                   )
///     (vertical regression matrix bin_m)
/// Since each matrix is square we can easily infer the bin index from the row index
/// First index of vertRegView is bin_index * number of levels + level index,
///     the second is number of levels associated with matrix column.
///
/// The interpolation weights are calculated for each grid point location
/// The second index relates to the bin index.
/// We ensure that across all bins for a grid point we sum to 1.
///
atlas::FieldSet createGpRegressionStats(const atlas::FunctionSpace & functionSpace,
                                        const atlas::FieldSet & fields,
                                        const oops::Variables & variables,
                                        const GpToHpCovarianceParameters & params) {
  // Get necessary parameters
  // path to covariance file with gp covariance parameters.
  const std::string covFileName(params.covariance_file_path);
  // number of latitudes that existed in the generation of the covariance file
  const std::size_t covGlobalNLats(static_cast<std::size_t>(params.covariance_nlat));
  // number of model levels
  const std::size_t modelLevels = variables["unbalanced_pressure_levels_minus_one"].getLevels();
  // geostrophic pressure vertical regression statistics are grouped
  // into overlapping bins based on latitude;
  // number of bins associated with the gP vertical regression
  const std::size_t gPBins(static_cast<std::size_t>(params.gp_regression_bins));

  oops::Log::info() <<
    "gp regression no of bins = " << gPBins << std::endl;
  oops::Log::info() <<
    "gp regression model levels = " << modelLevels << std::endl;

  atlas::FieldSet gpStatistics;

  // If gPBins is 0 - then we switch off the balanced pressure contribution.
  if (gPBins > 0) {
    gpStatistics.add(createGpRegressionMatrices(covFileName, gPBins, modelLevels));

    gpStatistics.add(createGpRegressionWeights(functionSpace, fields,
                                               covFileName, covGlobalNLats, gPBins));
  }

  return gpStatistics;
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
