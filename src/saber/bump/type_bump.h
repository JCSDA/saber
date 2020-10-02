/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_TYPE_BUMP_H_
#define SABER_BUMP_TYPE_BUMP_H_

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
                         atlas::functionspace::FunctionSpaceImpl *, atlas::field::FieldSetImpl *,
                         const eckit::Configuration &,
                         const eckit::Configuration &);
  void bump_add_member_f90(const int &, atlas::field::FieldSetImpl *, const int &, const int &);
  void bump_remove_member_f90(const int &, atlas::field::FieldSetImpl *, const int &, const int &);
  void bump_run_drivers_f90(const int &);
  void bump_apply_vbal_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_vbal_inv_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_vbal_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_vbal_inv_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_stddev_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_stddev_inv_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_apply_nicas_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_get_cv_size_f90(const int &, int &);
  void bump_apply_nicas_sqrt_f90(const int &, const double *, atlas::field::FieldSetImpl *);
  void bump_apply_nicas_sqrt_ad_f90(const int &, atlas::field::FieldSetImpl *, const double *);
  void bump_randomize_f90(const int &, atlas::field::FieldSetImpl *);
  void bump_get_parameter_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
  void bump_set_parameter_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
  void bump_dealloc_f90(const int &);
}
}  // namespace saber

#endif  // SABER_BUMP_TYPE_BUMP_H_
