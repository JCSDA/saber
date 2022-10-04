/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit {
  class LocalConfiguration;
}

namespace saber {
  namespace gsi {
    typedef int GridKey;
    extern "C" {
      void gsi_grid_create_f90(GridKey &, const eckit::Configuration &,
                              const eckit::mpi::Comm &);
      void gsi_grid_delete_f90(GridKey &);
      void gsi_grid_print_f90(const GridKey &);
      void gsi_grid_get_levels_f90(const GridKey &, int &);
      void gsi_grid_set_atlas_lonlat_f90(const GridKey &, atlas::field::FieldSetImpl *);
    }
  }  // namespace gsi
}  // namespace saber
