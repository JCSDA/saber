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

  // Change variables in the background to inner variables
  // TODO(AS): this code should also be in calibrationInverseMultiply
  oops::Variables neededVars = innerVars_;
  atlas::FieldSet xb_inner = xb;
  oops::Variables varsVaderPopulates = vader_.changeVar(xb_inner, neededVars);
  ASSERT(varsVaderPopulates == innerVars_);
  // Pass everything but outer variables to the vader TL/AD execution plan
  // (xb_outer should still have all inner variables)
  atlas::FieldSet xb_outer;
  for (const auto & fieldname : xb_inner.field_names()) {
     if (!outerVars_.has(fieldname) || innerVars_.has(fieldname)) {
       xb_outer.add(xb_inner[fieldname]);
     }
  }
  neededVars = outerVars_;
  varsVaderPopulates = vader_.changeVarTraj(xb_outer, neededVars);
  ASSERT(varsVaderPopulates == outerVars_);

  oops::Log::trace() << classname() << "::VaderBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Variables vars = outerVars_;
  vader_.changeVarTL(fset, vars);
  // copy only outer variables to the output fieldset (vader leaves both
  // output and input variables in the fieldset)
  atlas::FieldSet fset_out;
  for (const auto & fieldname : outerVars_.variables()) {
    fset_out.add(fset[fieldname]);
  }
  fset = fset_out;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Variables vars = outerVars_;
  vader_.changeVarAD(fset, vars);
  // copy only inner variables to the output fieldset (vader leaves both
  // output and input variables in the fieldset)
  atlas::FieldSet fset_out;
  for (const auto & fieldname : innerVars_.variables()) {
    fset_out.add(fset[fieldname]);
  }
  fset = fset_out;
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
  os << "Vader linear variable change from " << innerVars_
     << " to " << outerVars_;
}

// -----------------------------------------------------------------------------

}  // namespace saber
