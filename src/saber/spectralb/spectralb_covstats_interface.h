/*
 * (C) Crown Copyright 2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef SABER_SPECTRALB_SPECTRALB_COVSTATS_INTERFACE_H_
#define SABER_SPECTRALB_SPECTRALB_COVSTATS_INTERFACE_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace saber {
namespace spectralb {

extern "C" {

void covSpectralBins_f90(
  const eckit::Configuration &,
  const int &,
  const char *,
  const int &);

void covSpectralUMatrix_f90(
  const eckit::Configuration &,
  const int &,
  const char *,
  const int &,
  const int &,
  const int &,
  float &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
#endif  // SABER_SPECTRALB_SPECTRALB_COVSTATS_INTERFACE_H_
