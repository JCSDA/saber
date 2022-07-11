/* 
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#ifndef SABER_SPECTRALB_COVARIANCESTATISTICS_H_
#define SABER_SPECTRALB_COVARIANCESTATISTICS_H_

#include <Eigen/Core>
#include <Eigen/StdVector>

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

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"

#include "saber/spectralb/CovarianceStatisticsUtils.h"
#include "saber/spectralb/spectralbParameters.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// 'CovStat_ErrorCov': class for covariance statistics associated with the error covariance;
template<typename MODEL>
class CovStat_ErrorCov {
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef spectralbParameters<MODEL> Parameters_;

 public:
  static const std::string classname() {return "saber::CovStat_ErrorCov";}

  CovStat_ErrorCov(const Geometry_ & geom_,
                   const oops::Variables & vars,
                   const Parameters_ &);

  /// \details getSpectralUMatrix() gets the square root of the spectral vertical covariances for
  ///          each total wavenumber and active variable
  ///          spectral vertical covariance = UMatrix UMatrix^T
  const std::map<std::string, std::vector<Eigen::MatrixXf>> & getSpectralUMatrices()  const {
    return spectralUMatrices_;
  }

  /// \details getSpectralVerticalCovariances() gets the spectral vertical covariances for
  ///          each total wavenumber and active variable
  const std::map<std::string, std::vector<Eigen::MatrixXf>> &
      getSpectralVerticalCovariances()  const {
    return spectralVerticalCovariances_;
  }

  /// \details getSpectralVerticalCovariances() gets the spectral vertical covariances for
  ///          each total wavenumber and active variable
  const std::map<std::string, std::vector<Eigen::MatrixXf>> &
      getSpectralVerticalCorrelations()  const {
    return spectralVerticalCorrelations_;
  }

 private:
  // name of covariance file
  std::string covarianceFileName_;
  // number of model levels
  int modelLevels_;
  // number of spectral bins for each field
  std::vector<std::size_t> netCDFSpectralBins_;
  // square root of the spectral vertical covariances
  // with the number of spectral bins will be that of the cov file
  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralUMatrices_;
  // spectral vertical covariances
  // with the number of spectral bins to be that for the gaussian grid resolution
  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralVerticalCovariances_;
  // spectral standard deviations
  std::map<std::string, Eigen::VectorXd> spectralSD_;
  // spectral vertical correlations
  // with the number of spectral bins to be that for the gaussian grid resolution
  std::map<std::string, std::vector<Eigen::MatrixXf>> spectralVerticalCorrelations_;

  void print(std::ostream &) const;
};

}  // namespace spectralb
}  // namespace saber
// end of Header definition

// start of Header implementation
namespace saber {
namespace spectralb {

template<typename MODEL>
CovStat_ErrorCov<MODEL>::CovStat_ErrorCov(const Geometry_ & geom_,
                                          const oops::Variables & vars,
                                          const Parameters_ & params) :
  covarianceFileName_(params.covarianceFile),
  modelLevels_(geom_.variableSizes(vars)[0]),
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
  for (std::size_t lvl : geom_.variableSizes(vars)) {
    if (static_cast<int>(lvl) != modelLevels_) {
      throw eckit::UnexpectedState("spectral covariance block assumes all fields have "
                                   "same number of model levels");
    }
  }
}

template<typename MODEL>
void CovStat_ErrorCov<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print starting" << std::endl;

  os << std::endl << "  covstats print";
  oops::Log::trace() <<
    "Covariance Statistics (SABER spectral B, error covariance) print done" << std::endl;
}

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_COVARIANCESTATISTICS_H_
