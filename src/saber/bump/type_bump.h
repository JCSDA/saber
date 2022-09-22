/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace saber {
extern "C" {
  void bump_create_f90(int &, const eckit::mpi::Comm *,
                       const atlas::functionspace::FunctionSpaceImpl *,
                       const atlas::field::FieldSetImpl *,
                       const eckit::Configuration &,
                       const eckit::Configuration &,
                       const atlas::field::FieldSetImpl *);
  void bump_second_geometry_f90(int &,
                                const atlas::functionspace::FunctionSpaceImpl *,
                                const atlas::field::FieldSetImpl *);
  void bump_add_member_f90(const int &, const atlas::field::FieldSetImpl *,
                           const int &, const int &);
  void bump_update_vbal_cov_f90(const int &, const atlas::field::FieldSetImpl *,
                                const int &);
  void bump_update_var_f90(const int &, const atlas::field::FieldSetImpl *,
                           const int &);
  void bump_update_mom_f90(const int &, const atlas::field::FieldSetImpl *,
                           const int &, const int &);
  void bump_run_drivers_f90(const int &);
  void bump_apply_vbal_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_vbal_inv_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_vbal_ad_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_vbal_inv_ad_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_stddev_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_stddev_inv_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_apply_nicas_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_get_cv_size_f90(const int &, int &);
  void bump_randomize_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_psichi_to_uv_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_psichi_to_uv_ad_f90(const int &, const atlas::field::FieldSetImpl *);
  void bump_get_parameter_f90(const int &, const int &, const char *,
                              const int &, const int &, const atlas::field::FieldSetImpl *);
  void bump_set_ncmp_f90(const int &, const int &, const int &);
  void bump_set_parameter_f90(const int &, const int &, const char *,
                              const int &, const atlas::field::FieldSetImpl *);
  void bump_partial_dealloc_f90(const int &);
  void bump_dealloc_f90(const int &);
}
}  // namespace saber
