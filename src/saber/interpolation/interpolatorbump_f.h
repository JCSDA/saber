/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INTERPOLATORBUMP_F_H_
#define OOPS_GENERIC_INTERPOLATORBUMP_F_H_

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
  void bint_create_f90(int &, const eckit::mpi::Comm *,
                       const atlas::functionspace::FunctionSpaceImpl *,
                       const atlas::functionspace::FunctionSpaceImpl *,
                       const atlas::field::FieldSetImpl *,
                       const eckit::Configuration * const *);
  void bint_delete_f90(const int &);
  void bint_apply_f90(const int &, const atlas::field::FieldSetImpl *,
                      atlas::field::FieldSetImpl *);
  void bint_apply_ad_f90(const int &, const atlas::field::FieldSetImpl *,
                         atlas::field::FieldSetImpl *);
}
}  // namespace oops

#endif  // OOPS_GENERIC_INTERPOLATORBUMP_F_H_
