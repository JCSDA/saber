/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierGridToSpectral.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierGridToSpectral>
  makerBifourierGridToSpectral_("BifourierGridToSpectral");

// -----------------------------------------------------------------------------

BifourierGridToSpectral::BifourierGridToSpectral(const oops::GeometryData & outerGeometryData,
                                                 const oops::Variables & outerVars,
                                                 const eckit::Configuration & covarConfig,
                                                 const Parameters_ & params,
                                                 const oops::FieldSet3D & xb,
                                                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerVars_(outerVars),
    trans_(transStore_.retrieveTransform(outerGeometryData)),
    innerGeometryData_(trans_->geometryData())
{
  oops::Log::trace() << classname() << "::BifourierGridToSpectral starting" << std::endl;
  oops::Log::trace() << classname() << "::BifourierGridToSpectral done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierGridToSpectral::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Direct spectral transform
  atlas::FieldSet fsetTmp;
  trans_->gp2sp(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Copy FieldSet
  trans_->copyFieldSet(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierGridToSpectral::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Copy FieldSet
  atlas::FieldSet fsetTmp;
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Direct spectral transform, adjoint
  trans_->gp2spAdj(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierGridToSpectral::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  // Copy FieldSet
  atlas::FieldSet fsetTmp;
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp, innerVars_);

  // Remove variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Inverse spectral transform
  trans_->sp2gp(fsetTmp, fset.fieldSet(), innerVars_);

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierGridToSpectral::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
