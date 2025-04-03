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
  // In the following code, "inner" and "outer" are used to refer to the two sets of variables,
  // according to how close they are to the center of the linear operator chain. These are the
  // least ambiguous names.
  // However, when maintaining this code, sometimes regular humans may find it helpful to think of
  // them with other names (which might not always be completely accurate):
  // "inner" ~ "control" ~ "B-matrix" ~ "from" ~ "ingredients"
  // "outer" ~ "increment" ~ "analysis" ~ "background" ~ "to" ~ "products"

  oops::Log::trace() << classname() << "::VaderBlock starting" << std::endl;
  // Change variables in the background to inner variables
  // TODO(someone): perhaps this code will happen in the ErrorCovariance ctor?
  oops::Variables neededVars = outerVars_;
  atlas::FieldSet xb_outer = xb.fieldSet();
  oops::Variables ingredientVars = innerVars_;

  // We pass xb_outer to Vader to store in its trajectory, even though the trajectory "should"
  // contain inner variables. However, in initTLAD, Vader will attempt to produce any (inner)
  // trajectory variables that are required by its linear recipes, so there's no need for us
  // to try to create them in advance here.
  vader_.changeVarTraj(xb_outer, neededVars);
  const oops::Variables varsProduced = vader_.initTLAD(ingredientVars);
  neededVars -= varsProduced;
  ASSERT_MSG(neededVars.size() == 0, "Vader could not produce all outer variables");

  oops::Log::trace() << classname() << "::VaderBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

void VaderBlock::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  vader_.changeVarTL(fset.fieldSet());
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
  vader_.changeVarAD(fset.fieldSet());
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
