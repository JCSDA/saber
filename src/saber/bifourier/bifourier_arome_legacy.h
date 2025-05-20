/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

extern "C" {
  void bifourier_arome_legacy_vortopb_f90(const eckit::Configuration &,
                                          const int &,
                                          double[]);

  void bifourier_arome_legacy_balance_f90(const eckit::Configuration &,
                                          const int &,
                                          const int &,
                                          double[],
                                          double[],
                                          double[],
                                          double[],
                                          double[],
                                          double[]);

  void bifourier_arome_legacy_covariance_f90(const eckit::Configuration &,
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
