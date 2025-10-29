/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/Biperiodization.h"

#include <algorithm>

#include "atlas/grid.h"

#include "oops/util/FloatCompare.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<Biperiodization> makerBiperiodization_("Biperiodization");

// -----------------------------------------------------------------------------

Biperiodization::Biperiodization(const oops::GeometryData & outerGeometryData,
                                 const oops::Variables & outerVars,
                                 const eckit::Configuration & covarConfig,
                                 const Parameters_ & params,
                                 const oops::FieldSet3D & xb,
                                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    outerGeometryData_(outerGeometryData),
    comm_(outerGeometryData_.comm()),
    innerVars_(outerVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::Biperiodization starting" << std::endl;

  // Setup biperiodization implementation
  biper_.reset(new BiperiodizationImpl(outerGeometryData_, outerVars, params_.biperParams.value()));

  // Empty inner FieldSet
  atlas::FieldSet innerFset;

  if (!biper_->sameGrid()) {
    // Generate inner GeometryData
    innerGeometryData_.reset(new oops::GeometryData(biper_->innerFunctionSpace(), innerFset,
      outerGeometryData_.levelsAreTopDown(), comm_));
  }

  if (params_.read.value() != boost::none) {
    // Create input test fieldset
    inputTestFset_.reset(new oops::FieldSet3D(xb.validTime(), outerGeometryData_.comm()));
  }

  oops::Log::trace() << classname() << "::Biperiodization done" << std::endl;
}

// -----------------------------------------------------------------------------

const oops::GeometryData & Biperiodization::innerGeometryData() const {
  if (innerGeometryData_) {
    return *innerGeometryData_;
  } else {
    return outerGeometryData_;
  }
}

// -----------------------------------------------------------------------------

void Biperiodization::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Apply multiply
  biper_->multiply(fset.fieldSet());

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Biperiodization::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Apply multiply adjoint
  biper_->multiplyAD(fset.fieldSet());

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void Biperiodization::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  // Apply multiply left inverse
  biper_->leftInverseMultiply(fset.fieldSet());

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void Biperiodization::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Read input test file
  inputTestFset_->read(outerGeometryData_.functionSpace(),
                       innerVars_,
                       params_.read.value()->inputTestFile.value());

  // LeftInverseMultiply on input test file
  outputInnerTestFset_.reset(new oops::FieldSet3D(*inputTestFset_));
  outputInnerTestFset_->name() = "outputInnerTestFile";
  leftInverseMultiply(*outputInnerTestFset_);

  // Multiply on output inner test file
  outputOuterTestFset_.reset(new oops::FieldSet3D(*outputInnerTestFset_));
  outputOuterTestFset_->name() = "outputOuterTestFile";
  multiply(*outputOuterTestFset_);

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void Biperiodization::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  const auto & paramsRead = params_.read.value();
  if (paramsRead != boost::none) {
    // Write output inner test file
    outputInnerTestFset_->write(paramsRead->outputInnerTestFile.value());

    // Write output outer file
    outputOuterTestFset_->write(paramsRead->outputOuterTestFile.value());
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void Biperiodization::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
