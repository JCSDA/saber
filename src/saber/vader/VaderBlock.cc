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

#include "oops/base/FieldSet3D.h"

#include "saber/oops/Utilities.h"


namespace saber {

namespace {

oops::Variables createInnerVars(
    const oops::Variables & outerVars,
    const oops::Variables & innerActiveVarsWithoutMeta) {
  // TO DO: Ideally we should get model levels of innerVars from vader.
  // For now we assume that
  // All vertical levels are the same in outerVars and that innerVars
  // has the same number of model levels
  oops::Variables innerActiveVarsWithMeta(innerActiveVarsWithoutMeta);
  std::vector<int> modelLevels;
  for (const auto & var : outerVars) {
    modelLevels.push_back(var.getLevels());
  }

  if (outerVars.size() > 1) {
     ASSERT(std::equal(modelLevels.begin() + 1, modelLevels.end(),
                       modelLevels.begin()));
  }

  for (auto & var : innerActiveVarsWithMeta) {
    var.setLevels(modelLevels[0]);
  }

  return innerActiveVarsWithMeta;
}

}  // namespace


// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VaderBlock> makerVaderBlock_("vader variable change");

// -----------------------------------------------------------------------------

VaderBlock::VaderBlock(const oops::GeometryData & outerGeometryData,
                       const oops::Variables & outerVars,
                       const eckit::Configuration & outerBlockConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    outerVars_(outerVars),
    innerGeometryData_(outerGeometryData),
    innerVars_(createInnerVars(outerVars, params.innerVars)),
    vader_(params.vader, outerBlockConf.getSubConfiguration("vader"))
{
  oops::Log::trace() << classname() << "::VaderBlock starting" << std::endl;

  // Change variables in the background to inner variables
  // TODO(someone): perhaps this code will happen in the ErrorCovariance ctor?
  oops::Variables neededVars = innerVars_;
  atlas::FieldSet xb_inner = xb.fieldSet();

  oops::Variables varsVaderPopulates = vader_.changeVar(xb_inner, neededVars);
  ASSERT_MSG(varsVaderPopulates == innerVars_, "VADER can not populate all "
             "inner variables for SABER block.");
  // Pass only inner variables to the vader TL/AD execution plan
  // (xb_inner has both outer and inner variables after calling vader::changeVar)
  atlas::FieldSet xb_outer;
  for (const auto & innerVar : innerVars_) {
    xb_outer.add(xb_inner[innerVar.name()]);
  }
  // Set trajectory and create a vader plan for going from inner to outer
  // variables.
  neededVars = outerVars_;
  varsVaderPopulates = vader_.changeVarTraj(xb_outer, neededVars);
  ASSERT(varsVaderPopulates == outerVars_);

  oops::Log::trace() << classname() << "::VaderBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Variables vars = outerVars_;
  vader_.changeVarTL(fset.fieldSet(), vars);
  // copy only outer variables to the output fieldset (vader leaves both
  // output and input variables in the fieldset)
  atlas::FieldSet fset_out;
  for (const auto & outerVar : outerVars_) {
    fset_out.add(fset[outerVar.name()]);
  }
  fset.fieldSet() = fset_out;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Variables vars = outerVars_;
  vader_.changeVarAD(fset.fieldSet(), vars);
  // copy only inner variables to the output fieldset (vader leaves both
  // output and input variables in the fieldset)
  atlas::FieldSet fset_out;
  for (const auto & innerVar : innerVars_) {
    fset_out.add(fset[innerVar.name()]);
  }
  fset.fieldSet() = fset_out;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::print(std::ostream & os) const {
  os << "Vader linear variable change from " << innerVars_
     << " to " << outerVars_;
}

// -----------------------------------------------------------------------------

}  // namespace saber
