/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/interpolation/Interpolation.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/grid/Grid.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberOuterBlockMaker<Interpolation> makerInterpolation_("gsi interpolation to model grid");

// -------------------------------------------------------------------------------------------------

Interpolation::Interpolation(const oops::GeometryData & outerGeometryData,
                             const std::vector<size_t> & activeVariableSizes,
                             const oops::Variables & outerVars,
                             const Parameters_ & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const std::vector<atlas::FieldSet> & fsetVec)
  : innerVars_(outerVars)
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;
  util::Timer timer(classname(), "Interpolation");

  // Grid
  Grid grid(outerGeometryData.comm(), params.toConfiguration());

  // Inner geometry and variables
  innerGeometryData_.reset(new oops::GeometryData(grid.functionSpace(),
                                                  outerGeometryData.fieldSet(),
                                                  outerGeometryData.levelsAreTopDown(),
                                                  outerGeometryData.comm()));

  // Active variables
  std::vector<std::string> activeVars =
    params.activeVars.value().get_value_or(outerVars).variables();

  // Create the interpolator
  interpolator_.reset(new UnstructuredInterpolation(outerGeometryData.comm(),
                                                    params.toConfiguration(),
                                                    innerGeometryData_->functionSpace(),
                                                    outerGeometryData.functionSpace(),
                                                    activeVariableSizes,
                                                    activeVars));

  // Create the interpolator
  inverseInterpolator_.reset(new UnstructuredInterpolation(outerGeometryData.comm(),
                                                           params.toConfiguration(),
                                                           outerGeometryData.functionSpace(),
                                                           innerGeometryData_->functionSpace(),
                                                           activeVariableSizes,
                                                           activeVars));

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Interpolation::~Interpolation() {
  oops::Log::trace() << classname() << "::~Interpolation starting" << std::endl;
  util::Timer timer(classname(), "~Interpolation");
  oops::Log::trace() << classname() << "::~Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  interpolator_->apply(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  interpolator_->applyAD(fset);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "calibrationInverseMultiply");
  inverseInterpolator_->apply(fset);
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
