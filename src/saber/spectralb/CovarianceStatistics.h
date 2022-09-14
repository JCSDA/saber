/* 
 * (C) Crown Copyright 2017-2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#ifndef SABER_SPECTRALB_COVARIANCESTATISTICS_H_
#define SABER_SPECTRALB_COVARIANCESTATISTICS_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "atlas/field.h"

#include "saber/spectralb/spectralbParameters.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// 'CovStat_ErrorCov': class for covariance statistics associated with the error covariance;
class CovStat_ErrorCov {
  typedef spectralbParameters Parameters_;

 public:
  static const std::string classname() {return "saber::CovStat_ErrorCov";}

  CovStat_ErrorCov(const std::vector<size_t> &,
                   const oops::Variables &,
                   const Parameters_ &);

  /// \details getSpectralUMatrix() gets the square root of the spectral vertical covariances for
  ///          each total wavenumber and active variable
  ///          spectral vertical covariance = UMatrix UMatrix^T
  const atlas::FieldSet & getSpectralUMatrices()  const {
    return spectralUMatrices_;
  }

  /// \details getSpectralVerticalCovariances() gets the spectral vertical covariances for
  ///          each total wavenumber and active variable
  const atlas::FieldSet & getSpectralVerticalCovariances()  const {
    return spectralVerticalCovariances_;
  }

  /// \details getSpectralVerticalCovariances() gets the spectral vertical covariances for
  ///          each total wavenumber and active variable
  const atlas::FieldSet & getSpectralVerticalCorrelations()  const {
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
  atlas::FieldSet spectralUMatrices_;
  // spectral vertical covariances
  // with the number of spectral bins to be that for the gaussian grid resolution
  atlas::FieldSet spectralVerticalCovariances_;
  // spectral standard deviations with model level
  atlas::FieldSet spectralSD_;
  // spectral vertical correlations
  // with the number of spectral bins to be that for the gaussian grid resolution
  atlas::FieldSet spectralVerticalCorrelations_;

  void print(std::ostream &) const;
};

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_COVARIANCESTATISTICS_H_
