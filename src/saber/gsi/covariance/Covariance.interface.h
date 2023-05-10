/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit {
  class LocalConfiguration;
}

namespace saber {
  namespace gsi {
    typedef int CovarianceKey;
    extern "C" {
      void gsi_covariance_create_f90(CovarianceKey &, const eckit::mpi::Comm &,
                                     const eckit::Configuration &,
                                     const atlas::field::FieldSetImpl *,
                                     const atlas::field::FieldSetImpl *,
                                     const util::DateTime &);
      void gsi_covariance_delete_f90(CovarianceKey &);
      void gsi_covariance_randomize_f90(const CovarianceKey &,
                                        const atlas::field::FieldSetImpl *);
      void gsi_covariance_multiply_f90(const CovarianceKey &, const atlas::field::FieldSetImpl *);
    }
  }  // namespace gsi
}  // namespace saber
