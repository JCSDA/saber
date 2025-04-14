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
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::Biperiodization starting" << std::endl;

  // Setup biperiodization implementation
  biper_.reset(new BiperiodizationImpl(outerGeometryData, outerVars, params_.biperParams.value()));

  // Empty inner FieldSet
  atlas::FieldSet innerFset;

  // Generate inner GeometryData
  innerGeometryData_.reset(new oops::GeometryData(biper_->innerFunctionSpace(), innerFset,
    outerGeometryData.levelsAreTopDown(), comm_));

  oops::Log::trace() << classname() << "::Biperiodization done" << std::endl;
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

std::vector<std::pair<std::string, eckit::LocalConfiguration>> Biperiodization::getReadConfs()
  const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> pairs;
  const auto & paramsRead = params_.read.value();
  if (paramsRead != boost::none) {
    pairs.push_back(std::pair<std::string, eckit::LocalConfiguration>("inputTestFile",
      paramsRead->inputTestFile.value()));
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void Biperiodization::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  if (fsetVec.size() == 1) {
    // Copy input test file
    inputTestFset_.reset(new oops::FieldSet3D(fsetVec[0]));
  } else if (fsetVec.size() > 1) {
    // Unexpected
    throw eckit::Exception("wrong number of fields to read", Here());
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> Biperiodization::fieldsToWrite()
  const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> pairs;
  const auto & paramsRead = params_.read.value();
  if (paramsRead != boost::none) {
    pairs.push_back(std::pair<eckit::LocalConfiguration, oops::FieldSet3D>(
      paramsRead->outputInnerTestFile.value(), *outputInnerTestFset_));
    pairs.push_back(std::pair<eckit::LocalConfiguration, oops::FieldSet3D>(
      paramsRead->outputOuterTestFile.value(), *outputOuterTestFset_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void Biperiodization::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Check that the input test file has been read
  ASSERT(inputTestFset_);

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

void Biperiodization::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
