/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/gsi/interpolation/GSI_InterpolationImpl.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class InterpolationParameters : public InterpolationImplParameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, InterpolationImplParameters)
};

// -------------------------------------------------------------------------------------------------

class Interpolation : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::gsi::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const eckit::mpi::Comm &,
                const atlas::FunctionSpace &,
                const atlas::FieldSet &,
                const std::vector<size_t> &,
                const Parameters_ &,
                const atlas::FieldSet &,
                const atlas::FieldSet &,
                const std::vector<atlas::FieldSet> &);
  virtual ~Interpolation();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  // GSI Interpolation implementation
  InterpolationImpl interpolationImpl_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
