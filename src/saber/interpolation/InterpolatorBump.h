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

#include "oops/base/InterpolatorBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"

namespace saber {

// -----------------------------------------------------------------------------

/*! \brief Interface for Bump interpolation
 *
 */

class InterpolatorBump : public oops::InterpolatorBase,
                         private util::ObjectCounter<InterpolatorBump> {
 public:
  static const std::string classname() {return "saber::InterpolatorBump";}

  InterpolatorBump(const eckit::Configuration &,
                   const atlas::FunctionSpace &,
                   const atlas::FunctionSpace &,
                   const atlas::field::FieldSetImpl * = nullptr,
                   const eckit::mpi::Comm & = oops::mpi::world());

  ~InterpolatorBump();

  void apply(atlas::FieldSet const &, atlas::FieldSet &) override;
  void apply_ad(atlas::FieldSet const &, atlas::FieldSet &) override;

 private:
  int keyBumpInterpolator_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------


// gnu compilers want the Bump interpolator factory to be defined here.
static oops::InterpolatorMaker<InterpolatorBump> makerBumpInterpolator_("bump");

}  // namespace saber

#endif  // SABER_INTERPOLATION_INTERPOLATORBUMP_H_
