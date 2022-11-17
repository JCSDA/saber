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

#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"

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

  InterpolationImpl(const Comm_ &,
                    const eckit::Configuration &,
                    const PointCloud_ &,
                    const PointCloud_ &,
                    const std::vector<size_t> &,
                    const std::vector<std::string> &);
  ~InterpolationImpl();
  void multiply(atlas::FieldSet &);
  void multiplyAD(atlas::FieldSet &);

 private:
  void print(std::ostream &) const;
  // Interpolation object
  std::unique_ptr<UnstructuredInterpolation> interpolator_;
  // Fieldsets
  atlas::FieldSet innerFields_;
  atlas::FieldSet outerFields_;
};

// -------------------------------------------------------------------------------------------------

InterpolationImpl::InterpolationImpl(const Comm_ & comm,
                                     const eckit::Configuration & conf,
                                     const PointCloud_ & innerFuncSpace,
                                     const PointCloud_ & outerFuncSpace,
                                     const std::vector<size_t> & activeVariableSizes,
                                     const std::vector<std::string> & activeVars)
  : interpolator_(), innerFields_(atlas::FieldSet()), outerFields_(atlas::FieldSet())
{
  oops::Log::trace() << classname() << "::InterpolationImpl starting" << std::endl;
  util::Timer timer(classname(), "InterpolationImpl");

  // Create the interpolator
  interpolator_.reset(new UnstructuredInterpolation(conf,
                                                    innerFuncSpace,
                                                    outerFuncSpace,
                                                    activeVars,
                                                    comm));

  // Create the intermediate FieldSets to save repeated allocation
  for (size_t i = 0; i < activeVars.size(); ++i) {
    innerFields_.add(innerFuncSpace.createField<double>(
      atlas::option::name(activeVars[i]) | atlas::option::levels(activeVariableSizes[i])));
    outerFields_.add(outerFuncSpace.createField<double>(
      atlas::option::name(activeVars[i]) | atlas::option::levels(activeVariableSizes[i])));
  }

  oops::Log::trace() << classname() << "::InterpolationImpl done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

InterpolationImpl::~InterpolationImpl() {
  oops::Log::trace() << classname() << "::~InterpolationImpl starting" << std::endl;
  util::Timer timer(classname(), "~InterpolationImpl");
  oops::Log::trace() << classname() << "::~InterpolationImpl done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void InterpolationImpl::multiply(atlas::FieldSet & fset) {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Do the interpolation from GSI grid to model grid
  interpolator_->apply(fset, outerFields_);

  // Replace the saber (model) fields with the GSI fields
  fset = outerFields_;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void InterpolationImpl::multiplyAD(atlas::FieldSet & fset) {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");

  // Do the adjoint of interpolation from GSI grid to model grid
  interpolator_->apply_ad(fset, innerFields_);

  // Replace the saber (model) fields with the GSI fields
  fset = innerFields_;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void InterpolationImpl::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
