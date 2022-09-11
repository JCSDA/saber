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

using atlas::option::levels;
using atlas::option::name;

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class InterpolationImplParameters : public GridParameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationImplParameters, GridParameters)

 public:
  // Interpolation method
  oops::Parameter<std::string> interpMethod{"interpolation method", "barycent", this};
};

// -------------------------------------------------------------------------------------------------

class InterpolationImpl {
 public:
  static const std::string classname() {return "saber::gsi::InterpolationImpl";}

  typedef atlas::functionspace::PointCloud PointCloud_;
  typedef eckit::mpi::Comm Comm_;
  typedef InterpolationImplParameters Parameters_;

  InterpolationImpl(const Comm_ &, const PointCloud_ &, const Parameters_ &,
                    const std::vector<std::string>);
  ~InterpolationImpl();
  void multiply(atlas::FieldSet &) const;
  void multiplyAD(atlas::FieldSet &) const;

 private:
  void print(std::ostream &) const;
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

InterpolationImpl::InterpolationImpl(const Comm_ & comm, const PointCloud_ & modelGrid,
                                     const Parameters_ & params,
                                     const std::vector<std::string> variables)
  : interpolator_(), variables_(), grid_(comm, params), modGridFuncSpace_(modelGrid)
{
  oops::Log::trace() << classname() << "::InterpolationImpl starting" << std::endl;
  util::Timer timer(classname(), "InterpolationImpl");

  // Object wide copy of the variables
  variables_ = variables;

  // Get number of levels
  gsiLevels_ = grid_.levels();

  // Get functionspace for the GSI grid
  gsiGridFuncSpace_ = atlas::functionspace::PointCloud(grid_.functionSpace());

  // Create the interpolator
  interpolator_.reset(new UnstructuredInterpolation(params.toConfiguration(),
                                                    gsiGridFuncSpace_,
                                                    modGridFuncSpace_,
                                                    nullptr, comm));

  oops::Log::trace() << classname() << "::InterpolationImpl done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

InterpolationImpl::~InterpolationImpl() {
  oops::Log::trace() << classname() << "::~InterpolationImpl starting" << std::endl;
  util::Timer timer(classname(), "~InterpolationImpl");
  oops::Log::trace() << classname() << "::~InterpolationImpl done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void InterpolationImpl::multiply(atlas::FieldSet & fset) const {
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

void InterpolationImpl::multiplyAD(atlas::FieldSet & fset) const {
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

void InterpolationImpl::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
