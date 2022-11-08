/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/VaderBlock.h"

#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"


namespace saber {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VaderBlock> makerVaderBlock_("vader variable change");

// -----------------------------------------------------------------------------

VaderBlock::VaderBlock(const oops::GeometryData & outerGeometryData,
                       const std::vector<size_t> & activeVariableSizes,
                       const oops::Variables & outerVars,
                       const Parameters_ & params,
                       const atlas::FieldSet & xb,
                       const atlas::FieldSet & fg,
                       const std::vector<atlas::FieldSet> & fsetVec)
  : outerVars_(outerVars),
    innerGeometryData_(outerGeometryData), innerVars_(params.innerVars),
    vader_(params.vader)
{
  oops::Log::trace() << classname() << "::VaderBlock starting" << std::endl;

  oops::Variables neededVars = innerVars_;
  // note: vader_.changeVarTraj changes the first argument which is why copy
  // is needed. would a copy like this imply that xb also gets changed (no
  // deep copy in atlas?)
  atlas::FieldSet xb_copy = xb;
  oops::Variables varsVaderPopulates = vader_.changeVarTraj(xb_copy, neededVars);
  assert(varsVaderPopulates == innerVars_);

  oops::Log::trace() << classname() << "::VaderBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Variables vars = innerVars_;
  vader_.changeVarTL(fset, vars);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Variables vars = innerVars_;
  vader_.changeVarAD(fset, vars);
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::print(std::ostream & os) const {
  os << "Vader linear variable change from " << outerVars_
     << " to " << innerVars_;
}

// -----------------------------------------------------------------------------

}  // namespace saber
