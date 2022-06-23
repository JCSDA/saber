/*
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef SABER_SPECTRALB_COVARIANCESTATISTICSUTILS_H_
#define SABER_SPECTRALB_COVARIANCESTATISTICSUTILS_H_

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "saber/spectralb/spectralb_covstats_interface.h"
#include "saber/spectralb/spectralbParameters.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace saber {
namespace spectralb {

template<typename MODEL>
std::vector<std::size_t> getNetCDFSpectralBins(
                             const spectralbParameters<MODEL> & params) {
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

template<typename MODEL>
std::map<std::string, std::vector<Eigen::MatrixXf>>
createUMatrices(const oops::Variables & activeVars,
                const int modelLevels,
                const std::vector<std::size_t> & netCDFSpectralBins,
                const spectralbParameters<MODEL> & params)
{
  // we are currently not using globalNons_ here
  // however in the near future we may want to scale the spectral covariances
  // so that they can be run at lower
  // spectral resolutions.  That is the reason for the globalNLons argument.

  oops::Variables netCDFVars(params.umatrixNetCDFNames);

  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralUMatrices;

  std::vector<Eigen::MatrixXf> spectral;

  for (std::size_t i = 0; i < activeVars.size(); ++i) {
    std::string var = activeVars[i];
    std::string netCDFVar = netCDFVars[i];

    std::vector<Eigen::MatrixXf> spectralUMatrix;

    // vector size
    int sizeVec = modelLevels * modelLevels * static_cast<int>(netCDFSpectralBins[i]);

    std::vector<float> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);

    covSpectralUMatrix_f90(params.toConfiguration(),
                           static_cast<int>(netCDFVar.size()),
                           netCDFVar.c_str(),
                           static_cast<int>(modelLevels),
                           static_cast<int>(netCDFSpectralBins[i]),
                           sizeVec,
                           spectralUMatrix1D[0]);

    Eigen::MatrixXf vertUvData(modelLevels, modelLevels);

    std::size_t jn(0);
    for (std::size_t b = 0; b < netCDFSpectralBins[i]; ++b) {
      for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(modelLevels); ++j) {
        for (Eigen::Index k = 0; k < static_cast<Eigen::Index>(modelLevels); ++k, ++jn) {
          vertUvData(j, k) = spectralUMatrix1D[jn];
        }
      }
      spectralUMatrix.push_back(vertUvData);
    }

    spectralUMatrices[var] = spectralUMatrix;
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
template<typename MODEL>
std::map<std::string, std::vector<Eigen::MatrixXf>>
createSpectralCovariances(
    const oops::Variables & activeVars,
    const int modelLevels,
    const std::vector<std::size_t> & netCDFSpectralBins,
    const std::map<std::string, std::vector<Eigen::MatrixXf>> & spectralUMatrices,
    const spectralbParameters<MODEL> & params)
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

  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralVerticalCovariances;

  std::vector<Eigen::MatrixXf> spectralVerticalCovariance;

  Eigen::MatrixXf vertCovData(modelLevels, modelLevels);

  for (std::size_t i = 0; i < activeVars.size(); ++i) {
    std::string var = activeVars[i];
    std::vector<Eigen::MatrixXf> spectralUMatrix =
      spectralUMatrices.at(var);

    std::vector<Eigen::MatrixXf> spectralVerticalCovariance;

    Eigen::MatrixXf vertCovData(modelLevels, modelLevels);

    for (int b = 0; b < nBins; ++b) {
      vertCovData = spectralUMatrix[b] *
                    spectralUMatrix[b].transpose();
      vertCovData = vertCovData * nBins / (netCDFSpectralBins[i]);
      spectralVerticalCovariance.push_back(vertCovData);
    }

    spectralVerticalCovariances[var] = spectralVerticalCovariance;
  }

  return spectralVerticalCovariances;
}

std::map<std::string, Eigen::VectorXd> createSpectralSD(
    const oops::Variables & activeVars,
    const int modelLevels,
    const std::map<std::string, std::vector<Eigen::MatrixXf>> & spectralVerticalCovariances)
{
  std::map<std::string, Eigen::VectorXd> spectralSDs;
  for (std::string var : activeVars.variables()) {
    std::vector<Eigen::MatrixXf> spectralVertCov =
      spectralVerticalCovariances.at(var);

    Eigen::VectorXd spectralSD(modelLevels);
    for (int k = 0; k < modelLevels; ++k) {
      spectralSD(k) = 0.0;
      for (std::size_t b = 0; b < spectralVertCov.size(); ++b) {
        spectralSD(k) += spectralVertCov[b](k, k);
      }
      spectralSD(k) = std::sqrt(spectralSD(k));
    }
    spectralSDs[var] = spectralSD;
  }

  return spectralSDs;
}

std::map<std::string, std::vector<Eigen::MatrixXf>> createSpectralCorrelations(
    const oops::Variables & activeVars,
    const int modelLevels,
    const std::map<std::string, std::vector<Eigen::MatrixXf>> & spectralVerticalCovariances,
    const std::map<std::string, Eigen::VectorXd> & spectralSDs)
{
  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralCorrelations;
  for (std::string var : activeVars.variables()) {
    std::vector<Eigen::MatrixXf> spectralVertCov =
      spectralVerticalCovariances.at(var);
    Eigen::VectorXd spectralSD = spectralSDs.at(var);

    std::vector<Eigen::MatrixXf> spectralCorrelation;

    Eigen::MatrixXf correlationScaling(modelLevels, modelLevels);
    for (int k1 = 0; k1 < modelLevels; ++k1) {
      for (int k2 = 0; k2 < modelLevels; ++k2) {
        correlationScaling(k1, k2) =
          static_cast<double>(spectralVertCov.size())/(spectralSD(k1) * spectralSD(k2));
      }
    }

    Eigen::MatrixXf correlation(modelLevels, modelLevels);

    for (std::size_t bin = 0; bin < spectralVertCov.size(); ++bin) {
      correlation = spectralVertCov[bin].cwiseProduct(correlationScaling);
      spectralCorrelation.push_back(correlation);
    }

    spectralCorrelations[var] = spectralCorrelation;
  }

  return spectralCorrelations;
}


}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_COVARIANCESTATISTICSUTILS_H_
