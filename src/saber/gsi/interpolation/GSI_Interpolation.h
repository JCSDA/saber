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
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/grid/GSI_Grid.h"
#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/oops/SaberBlockBase.h"

using atlas::option::levels;
using atlas::option::name;

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class InterpolationParameters : public GridParameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, GridParameters)

 public:
  // Interpolation method
  oops::Parameter<std::string> interpMethod{"interpolation method", "barycent", this};
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class Interpolation : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;

 public:
  static const std::string classname() {return "saber::gsi::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const Geometry_ &, const Parameters_ &, const State_ &, const State_ &);
  virtual ~Interpolation();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  // Pointer to Params
  std::unique_ptr<Parameters_> params_;
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

template<typename MODEL>
Interpolation<MODEL>::Interpolation(const Geometry_ & geom, const Parameters_ & params,
                                    const State_ & xb, const State_ & fg)
  : SaberBlockBase<MODEL>(params), params_(new Parameters_(params)),
    interpolator_(), variables_(), grid_(geom.getComm(), params)
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;
  util::Timer timer(classname(), "Interpolation");

  // Assert that there is no variable change in this block
  ASSERT(params.inputVars.value() == params.outputVars.value());
  variables_ = params.inputVars.value().variables();

  // Get number of levels
  gsiLevels_ = grid_.levels();

  // Get functionspace for the GSI grid
  gsiGridFuncSpace_ = atlas::functionspace::PointCloud(grid_.functionSpace());

  // Function space for the model grid
  modGridFuncSpace_ = atlas::functionspace::PointCloud(geom.functionSpace());

  // Create the interpolator
  interpolator_.reset(new UnstructuredInterpolation(params.toConfiguration(),
                                                    gsiGridFuncSpace_,
                                                    geom.functionSpace().get(),
                                                    nullptr, geom.getComm()));

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
Interpolation<MODEL>::~Interpolation() {
  oops::Log::trace() << classname() << "::~Interpolation starting" << std::endl;
  util::Timer timer(classname(), "~Interpolation");
  oops::Log::trace() << classname() << "::~Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::randomize(atlas::FieldSet & fset) const {
  ABORT(classname() + "randomize: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Create empty Model fieldset
  atlas::FieldSet modFields = atlas::FieldSet();

  // Loop over saber (gsi) fields and create corresponding model fields
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the model field and add to Fieldset
      modFields.add(modGridFuncSpace_.createField<double>(fieldName | levels(sabField.levels())));
  }

  // Do the interpolation from GSI grid to model grid
  interpolator_->apply(fset, modFields);

  // Replace the saber (model) fields with the GSI fields
  fset = modFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::inverseMultiply(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiply: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");

  // Create empty GSI fieldset
  atlas::FieldSet gsiFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding GSI fields
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Assert level match
      if (sabField.levels() != 0) {
        ASSERT(sabField.levels() == gsiLevels_);
      }

      // Create the field and add to Fieldset
      gsiFields.add(gsiGridFuncSpace_.createField<double>(fieldName | levels(sabField.levels())));
  }

  // Do the adjoint of interpolation from GSI grid to model grid
  interpolator_->apply_ad(fset, gsiFields);

  // Replace the saber (model) fields with the GSI fields
  fset = gsiFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::inverseMultiplyAD(atlas::FieldSet & fset) const {
  ABORT(classname() + "inverseMultiplyAD: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Interpolation<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
