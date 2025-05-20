/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierSplitTPs.h"

#include "atlas/field.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierSplitTPs> makerBifourierSplitTPs_("BifourierSplitTPs");

// -----------------------------------------------------------------------------

BifourierSplitTPs::BifourierSplitTPs(const oops::GeometryData & outerGeometryData,
                                     const oops::Variables & outerVars,
                                     const eckit::Configuration & covarConfig,
                                     const Parameters_ & params,
                                     const oops::FieldSet3D & xb,
                                     const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::BifourierSplitTPs starting" << std::endl;

  // Inner variables
  if (!params_.backward.value()) {
    // Add TPs to inner variables and remove T and Ps
    nz_ = innerVars_["air_temperature"].getLevels();
    innerVars_.push_back("air_temperature_and_log_of_air_pressure_at_surface");
    innerVars_["air_temperature_and_log_of_air_pressure_at_surface"].setLevels(nz_+1);
    innerVars_ -= innerVars_["air_temperature"];
    innerVars_ -= innerVars_["log_of_air_pressure_at_surface"];
  } else {
    // Add T and Ps to inner variables and remove TPs
    nz_ = innerVars_["air_temperature_and_log_of_air_pressure_at_surface"].getLevels()-1;
    innerVars_.push_back("air_temperature");
    innerVars_["air_temperature"].setLevels(nz_);
    innerVars_.push_back("log_of_air_pressure_at_surface");
    innerVars_["log_of_air_pressure_at_surface"].setLevels(1);
    innerVars_ -= innerVars_["air_temperature_and_log_of_air_pressure_at_surface"];
  }

  // Retrieve spectral transform
  trans_ = transStore_.retrieveTransform(outerGeometryData);

  oops::Log::trace() << classname() << "::BifourierSplitTPs done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Forward application
    forward(fset);
  } else {
    // Backward application
    backward(fset);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (!params_.backward.value()) {
    // Backward application
    backward(fset);
  } else {
    // Forward application
    forward(fset);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Backward application
    backward(fset);
  } else {
    // Forward application
    forward(fset);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::forward(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::forward starting" << std::endl;

  // Get inner field
  const auto tPsField = fset["air_temperature_and_log_of_air_pressure_at_surface"];

  // Create outer fields
  atlas::Field tField = trans_->spFspace()->createField<double>(
    atlas::option::name("air_temperature") | atlas::option::levels(nz_));
  atlas::Field psField = trans_->spFspace()->createField<double>(
    atlas::option::name("log_of_air_pressure_at_surface") | atlas::option::levels(1));

  // Get fields views
  const auto tPsView = make_view<double, 2>(tPsField);
  auto tView = make_view<double, 2>(tField);
  auto psView = make_view<double, 2>(psField);

  // Copy data
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      tView(js, jz) = tPsView(js, jz);
    }
    psView(js, 0) = tPsView(js, nz_);
  }

  // Remove inner field
  util::removeFieldsFromFieldSet(fset.fieldSet(),
    {"air_temperature_and_log_of_air_pressure_at_surface"});

  // Add outer fields
  fset.add(tField);
  fset.add(psField);

  oops::Log::trace() << classname() << "::forward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::backward(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::backward starting" << std::endl;

  // Get outer fields
  const auto tField = fset["air_temperature"];
  const auto psField = fset["log_of_air_pressure_at_surface"];

  // Create inner field
  atlas::Field tPsField = trans_->spFspace()->createField<double>(
    atlas::option::name("air_temperature_and_log_of_air_pressure_at_surface") |
    atlas::option::levels(nz_+1));

  // Get fields views
  const auto tView = make_view<double, 2>(tField);
  const auto psView = make_view<double, 2>(psField);
  auto tPsView = make_view<double, 2>(tPsField);

  // Copy data
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      tPsView(js, jz) = tView(js, jz);
    }
    tPsView(js, nz_) = psView(js, 0);
  }

  // Remove outer fields
  util::removeFieldsFromFieldSet(fset.fieldSet(), {"air_temperature",
    "log_of_air_pressure_at_surface"});

  // Add inner field
  fset.add(tPsField);

  oops::Log::trace() << classname() << "::backward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSplitTPs::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
