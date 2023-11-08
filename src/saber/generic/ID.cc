/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/ID.h"

#include <vector>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<IDCentral> makerIDCentral_("ID");

// -----------------------------------------------------------------------------

IDCentral::IDCentral(const oops::GeometryData & geometryData,
                     const oops::Variables & activeVars,
                     const eckit::Configuration & covarConf,
                     const Parameters_ & params,
                     const oops::FieldSet3D & xb,
                     const oops::FieldSet3D & fg) :
    SaberCentralBlockBase(params),
    geometryData_(geometryData),
    activeVars_(activeVars),
    ctlVecSize_(0)
{
  oops::Log::trace() << classname() << "::IDCentral starting" << std::endl;

  // Compute total number of levels
  size_t nlev = 0;
  for (const std::string & var : activeVars.variables()) {
    nlev += activeVars.getLevels(var);
  }

  // Compute control vector size
  ctlVecSize_ = nlev*geometryData_.functionSpace().size();

  oops::Log::trace() << classname() << "::IDCentral done" << std::endl;
}

// -----------------------------------------------------------------------------

IDCentral::~IDCentral() {
  oops::Log::trace() << classname() << "::~IDCentral starting" << std::endl;
  oops::Log::trace() << classname() << "::~IDCentral done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDCentral::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  for (const auto & var : activeVars_.variables()) {
      ASSERT(fset.has(var));
  }

  // Overwrite input fieldSet with random numbers
  const atlas::FieldSet newFieldSet = util::createRandomFieldSet(geometryData_.comm(),
                                                                 geometryData_.functionSpace(),
                                                                 activeVars_);

  for (const auto & var : activeVars_.variables()) {
    fset[var] = newFieldSet[var];
  }

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDCentral::multiplySqrt(const atlas::Field & cv,
                             atlas::FieldSet & fset,
                             const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Copy from control vector to fieldset
  size_t index = offset;
  const auto cvView = atlas::array::make_view<double, 1>(cv);
  for (auto & field : fset) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        view(jnode, jlevel) = cvView(index);
        ++index;
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDCentral::multiplySqrtAD(const atlas::FieldSet & fset,
                               atlas::Field & cv,
                               const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;

  // Copy from fieldset to control vector
  size_t index = offset;
  auto cvView = atlas::array::make_view<double, 1>(cv);
  for (const auto & field : fset) {
    const auto view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t jnode = 0; jnode < field.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        cvView(index) = view(jnode, jlevel);
        ++index;
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDCentral::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDCentral::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<IDOuter> makerIDOuter_("ID");

// -----------------------------------------------------------------------------

IDOuter::IDOuter(const oops::GeometryData & outerGeometryData,
                 const oops::Variables & outerVars,
                 const eckit::Configuration & covarConfig,
                 const Parameters_ & params,
                 const oops::FieldSet3D & xb,
                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars)
{
  oops::Log::trace() << classname() << "::IDOuter starting" << std::endl;
  oops::Log::trace() << classname() << "::IDOuter done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDOuter::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDOuter::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDOuter::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void IDOuter::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
