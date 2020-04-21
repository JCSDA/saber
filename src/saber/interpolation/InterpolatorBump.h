/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INTERPOLATORBUMP_H_
#define OOPS_GENERIC_INTERPOLATORBUMP_H_

#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/generic/Interpolator.h"
#include "oops/parallel/mpi/mpi.h"

namespace oops {

// -----------------------------------------------------------------------------

/*! \brief Interface for Bump interpolation
 * 
 */

class InterpolatorBump : public Interpolator,
                         private util::ObjectCounter<InterpolatorBump> {
 public:
  static const std::string classname() {return "oops::InterpolatorBump";}

  InterpolatorBump(const eckit::Configuration &,
                   const atlas::FunctionSpace &,
                   const atlas::FunctionSpace &,
                   const atlas::field::FieldSetImpl * = nullptr,
                   const eckit::mpi::Comm & = oops::mpi::comm());

  ~InterpolatorBump();

  void apply(atlas::FieldSet const &, atlas::FieldSet &);
  void apply_ad(atlas::FieldSet const &, atlas::FieldSet &);

 private:
  int keyBumpInterpolator_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_INTERPOLATORBUMP_H_
