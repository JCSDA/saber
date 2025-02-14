/*
 * (C) Crown Copyright 2022-2025 Met Office
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

extern "C" {

void covRegressionMatrices_f90(
  const int &,
  const char *,
  const int &,
  const int &,
  const int &,
  float &);

void covRegressionWeights_f90(
  const int &,
  const char *,
  const int &,
  const int &,
  const int &,
  int &,
  int &,
  float &,
  float &);

void covMuStats_f90(
  const int &,
  const char *,
  const int &,
  const char *,
  const int &,
  const int &,
  const int &,
  float &);

void oldMuStats_f90(
  const int &,
  const char *,
  const int &,
  const char *,
  const int &,
  const int &,
  double &);

void oldCovMIOStats_f90(
  const int &,
  const char *,
  const int &,
  const char *,
  const int &,
  const int &,
  double &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace saber

