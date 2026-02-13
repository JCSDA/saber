/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierID.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<BifourierID> makerBifourierID_("BifourierID");

// -----------------------------------------------------------------------------

BifourierID::BifourierID(const oops::GeometryData & geometryData,
                         const oops::Variables & centralVars,
                         const eckit::Configuration & covarConf,
                         const Parameters_ & params,
                         const oops::FieldSet3D & xb,
                         const oops::FieldSet3D & fg) :
    SaberCentralBlockBase(params, xb.validTime(), geometryData, centralVars),
    comm_(geometryData.comm()),
    trans_(transStore_.retrieveTransform(geometryData))
{
  oops::Log::trace() << classname() << "::BifourierID starting" << std::endl;
  oops::Log::trace() << classname() << "::BifourierID done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierID::~BifourierID() {
  oops::Log::trace() << classname() << "::~BifourierID starting" << std::endl;
  oops::Log::trace() << classname() << "::~BifourierID done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierID::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Create random spectral vector
  trans_->createRandomFieldSet(fset.fieldSet(), centralVars());

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierID::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Define control vector
  atlas::Field cv("cv", make_datatype<double>(), make_shape(ctlVecSize()));

  // No offset
  const size_t offset = 0;

  // Square-root adjoint multiply
  multiplySqrtAD(fset, cv, offset);

  // Square-root multiply
  multiplySqrt(cv, fset, offset);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierID::multiplySqrt(const atlas::Field & cv,
                               oops::FieldSet3D & fset,
                               const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Convert control vector to spectral FieldSet
  fset.clear();
  trans_->cv2fset(cv, fset.fieldSet(), centralVars(), offset);

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierID::multiplySqrtAD(const oops::FieldSet3D & fset,
                                 atlas::Field & cv,
                                 const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;

  // Convert spectral FieldSet to control vector
  trans_->fset2cv(fset.fieldSet(), cv, centralVars(), offset);

  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierID::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

