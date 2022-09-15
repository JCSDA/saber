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

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/gsi/grid/GSI_Grid.h"
#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class InterpolationImplParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(InterpolationImplParameters, SaberOuterBlockParametersBase)

 public:
  // Grid
  oops::RequiredParameter<GridParameters> grid{"grid", this};

  // Interpolation method
  oops::Parameter<std::string> interpMethod{"interpolation method", "barycent", this};
};

// -------------------------------------------------------------------------------------------------

class Interpolation : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::gsi::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const eckit::mpi::Comm &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const eckit::Configuration &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);
  virtual ~Interpolation();

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;

  // Interpolation object
  std::unique_ptr<UnstructuredInterpolation> interpolator_;
  // Function spaces
  atlas::FunctionSpace gsiGridFuncSpace_;
  atlas::FunctionSpace modGridFuncSpace_;
  // Variables
  std::vector<std::string> variables_;
  // Expected number of levels in GSI grid
  int gsiLevels_;
  // Grid
  Grid grid_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
