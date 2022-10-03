/* 
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "saber/spectralb/CovarianceStatistics.h"

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "atlas/array.h"
#include "atlas/field.h"

#include "saber/spectralb/CovarianceStatisticsUtils.h"
#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

CovStat_ErrorCov::CovStat_ErrorCov(const std::vector<size_t> & variableSizes,
                                   const oops::Variables & vars,
                                   const Parameters_ & params) :
  covarianceFileName_(params.covarianceFile),
  modelLevels_(variableSizes[0]),
  netCDFSpectralBins_(getNetCDFSpectralBins(params)),
  spectralUMatrices_(createUMatrices(vars, modelLevels_,
                                     netCDFSpectralBins_, params)),
  spectralVerticalCovariances_(createSpectralCovariances(
                               vars, modelLevels_, netCDFSpectralBins_,
                               spectralUMatrices_, params)),
  spectralSD_(createSpectralSD(vars, modelLevels_,
                               spectralVerticalCovariances_)),
  spectralVerticalCorrelations_(createSpectralCorrelations(
                                vars, modelLevels_, spectralVerticalCovariances_,
                                spectralSD_))
{
  for (std::size_t lvl : variableSizes) {
    if (static_cast<int>(lvl) != modelLevels_) {
      throw eckit::UnexpectedState("spectral covariance block assumes all fields have "
                                   "same number of model levels");
    }
  }
}

void CovStat_ErrorCov::print(std::ostream & os) const {
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print starting" << std::endl;

  os << std::endl << "  covstats print";
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print done" << std::endl;
}

}  // namespace spectralb
}  // namespace saber
