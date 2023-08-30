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
                     const oops::FieldSet3D & fg,
                     const size_t & timeRank) :
    SaberCentralBlockBase(params),
    geometryData_(geometryData),
    activeVars_(activeVars),
    timeRank_(timeRank)
{
  oops::Log::trace() << classname() << "::IDCentral starting" << std::endl;
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
                                                                 activeVars_,
                                                                 timeRank_);

  for (const auto & var : activeVars_.variables()) {
    fset[var] = newFieldSet[var];
  }

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
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
