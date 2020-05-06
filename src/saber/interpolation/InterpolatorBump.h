/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_INTERPOLATION_INTERPOLATORBUMP_H_
#define SABER_INTERPOLATION_INTERPOLATORBUMP_H_

#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/generic/Interpolator.h"
#include "oops/parallel/mpi/mpi.h"

namespace saber {

// -----------------------------------------------------------------------------

/*! \brief Interface for Bump interpolation
 *
 */

class InterpolatorBump : public oops::Interpolator,
                         private util::ObjectCounter<InterpolatorBump> {
 public:
  std::string classname() const override {return "saber::InterpolatorBump";}

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

}  // namespace saber

#endif  // SABER_INTERPOLATION_INTERPOLATORBUMP_H_
