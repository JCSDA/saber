/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include "eckit/config/Configuration.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

extern "C" {
  void bifourier_arome_legacy_read_balance_f90(const eckit::Configuration &,
                                               const int &,
                                               const int &,
                                               double[],
                                               double[],
                                               double[],
                                               double[],
                                               double[],
                                               double[],
                                               const int &,
                                               double[]);

  void bifourier_arome_legacy_write_balance_f90(const eckit::Configuration &,
                                                const int &,
                                                const int &,
                                                double[],
                                                double[],
                                                double[],
                                                double[],
                                                double[],
                                                double[],
                                                const int &,
                                                double[]);

  void bifourier_arome_legacy_read_covariance_f90(const eckit::Configuration &,
                                                  const int &,
                                                  const int &,
                                                  double[],
                                                  double[],
                                                  double[],
                                                  double[]);

  void bifourier_arome_legacy_write_covariance_f90(const eckit::Configuration &,
                                                   const int &,
                                                   const int &,
                                                   double[],
                                                   double[],
                                                   double[],
                                                   double[]);
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
