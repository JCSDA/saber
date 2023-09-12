/*
 * (C) Crown Copyright 2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace saber {
namespace spectralb {

extern "C" {

void covSpectralBinsLevels_f90(
  const eckit::Configuration &,
  const int &,
  const char *,
  int &,
  int &);

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
  float &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
