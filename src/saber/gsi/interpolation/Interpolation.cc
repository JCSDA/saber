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

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/gsi/grid/Grid.h"

#include "saber/oops/Utilities.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberOuterBlockMaker<Interpolation> makerInterpolation_("gsi interpolation to model grid");

// -------------------------------------------------------------------------------------------------

Interpolation::Interpolation(const oops::GeometryData & outerGeometryData,
                             const oops::Variables & outerVars,
                             const eckit::Configuration & covarConf,
                             const Parameters_ & params,
                             const oops::FieldSet3D & xb,
                             const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()), innerVars_(outerVars)
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
  const oops::Variables activeVars = getActiveVars(params, outerVars);
  std::vector<size_t> activeVariableSizes;
  for (const std::string & var : activeVars.variables()) {
    activeVariableSizes.push_back(activeVars.getLevels(var));
  }

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

void Interpolation::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  interpolator_->apply(fset.fieldSet());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  interpolator_->applyAD(fset.fieldSet());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "leftInverseMultiply");
  inverseInterpolator_->apply(fset.fieldSet());
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
