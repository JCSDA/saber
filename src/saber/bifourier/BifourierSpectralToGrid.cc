/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierSpectralToGrid.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierSpectralToGrid>
  makerBifourierSpectralToGrid_("BifourierSpectralToGrid");

// -----------------------------------------------------------------------------

BifourierSpectralToGrid::BifourierSpectralToGrid(const oops::GeometryData & outerGeometryData,
                                                 const oops::Variables & outerVars,
                                                 const eckit::Configuration & covarConfig,
                                                 const Parameters_ & params,
                                                 const oops::FieldSet3D & xb,
                                                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerVars_(outerVars),
    params_(params),
    trans_(transStore_.setupTransform(outerGeometryData, innerVars_, params_.transform.value()))
{
  oops::Log::trace() << classname() << "::BifourierSpectralToGrid starting" << std::endl;

  // Create inner GeometryData
  innerGeometryData_ = std::make_unique<oops::GeometryData>(trans_->spFspace(),
    outerGeometryData.fieldSet(), outerGeometryData.levelsAreTopDown(), outerGeometryData.comm(),
    false);

  oops::Log::trace() << classname() << "::BifourierSpectralToGrid done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralToGrid::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Copy FieldSet
  atlas::FieldSet fsetTmp;
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Inverse spectral transform
  trans_->sp2gp(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralToGrid::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Temporary fieldset
  atlas::FieldSet fsetTmp;

  // Inverse spectral transform, adjoint
  trans_->sp2gpAdj(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Copy FieldSet
  trans_->copyFieldSet(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralToGrid::inverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;

  // Direct spectral transform
  atlas::FieldSet fsetTmp;
  trans_->gp2sp(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Copy FieldSet
  trans_->copyFieldSet(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralToGrid::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
