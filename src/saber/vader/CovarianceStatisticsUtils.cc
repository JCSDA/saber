/*
 * (C) Crown Copyright 2022 Met Office
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
                                       const atlas::FieldSet & extraFields,
                                       const std::string & covFileName,
                                       const std::size_t covGlobalNLats,
                                       const std::size_t gpBins) {
  // get interpolation weights from file
  atlas::idx_t horizPts = extraFields[0].shape(0);

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
  auto interWgtFld = atlas::Field(std::string("interpolation_weights"),
    atlas::array::make_datatype<double>(),
    atlas::array::make_shape(horizPts, gpBins));

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

void populateInterpMuStats(atlas::FieldSet & augmentedStateFieldSet,
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
    return (i >= muBins - 1 ? muBins - 2 : i);
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
    for (int jl = 0; jl < augmentedStateFieldSet["rht"].levels(); ++jl) {
      double normField = normalisedField(static_cast<double>(RHtView(jn, jl)));
      int indx = index(normField, muBins);
      double w = weight(normField, indx);
      muView(jn, jl) = static_cast<double>(interp(w, covFldView(jl, indx),
                                                     covFldView(jl, indx + 1)));
    }
  }
  if (varName.compare("muA") == 0) {
    for (atlas::idx_t jn = 0; jn < augmentedStateFieldSet[varName].shape(0); ++jn) {
      for (int jl = 0; jl < augmentedStateFieldSet[varName].levels(); ++jl) {
        // note -  the actual mu inverse field in VAR involves an additional step
        //         that involves interpolation of a "mu_table".
        //         What is here is a simplification.
        muView(jn, jl) = 1.0 / muView(jn, jl);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
