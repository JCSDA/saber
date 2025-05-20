/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierVorDivToRedWind.h"

#include "atlas/field.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierVorDivToRedWind>
  makerBifourierVorDivToRedWind_("BifourierVorDivToRedWind");

// -----------------------------------------------------------------------------

BifourierVorDivToRedWind::BifourierVorDivToRedWind(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::BifourierVorDivToRedWind starting" << std::endl;

  // Inner variables
  if (!params.backward.value()) {
    // Add vor/div to inner variables and remove u/v
    nz_ = innerVars_["reduced_x_wind"].getLevels();
    innerVars_.push_back("air_upward_absolute_vorticity");
    innerVars_["air_upward_absolute_vorticity"].setLevels(nz_);
    innerVars_ -= innerVars_["reduced_x_wind"];
    innerVars_.push_back("air_horizontal_divergence");
    innerVars_["air_horizontal_divergence"].setLevels(nz_);
    innerVars_ -= innerVars_["reduced_y_wind"];
  } else {
    // Add u/v to inner variables and remove vor/div
    nz_ = innerVars_["air_upward_absolute_vorticity"].getLevels();
    innerVars_.push_back("reduced_x_wind");
    innerVars_["reduced_x_wind"].setLevels(nz_);
    innerVars_ -= innerVars_["air_upward_absolute_vorticity"];
    innerVars_.push_back("reduced_y_wind");
    innerVars_["reduced_y_wind"].setLevels(nz_);
    innerVars_ -= innerVars_["air_horizontal_divergence"];
  }

  // Retrieve spectral transform
  trans_ = transStore_.retrieveTransform(outerGeometryData);

  // Get grid function space
  const atlas::functionspace::StructuredColumns fs(trans_->geometryData().functionSpace());

  // Get js for (jk, jl) = (0, 0)
  jsZero_ = -1;
  for (size_t js = 0; js < trans_->ns(); ++js) {
    if ((trans_->jk(js) == 0) && (trans_->jl(js) == 0)) {
      jsZero_ = js;
    }
  }

  oops::Log::trace() << classname() << "::BifourierVorDivToRedWind done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::multiply(oops::FieldSet3D & fset) const {
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

void BifourierVorDivToRedWind::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (!params_.backward.value()) {
    // Forward application, adjoint
    forwardAD(fset);
  } else {
    // Backward application, adjoint
    backwardAD(fset);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::leftInverseMultiply(oops::FieldSet3D & fset) const {
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

void BifourierVorDivToRedWind::forward(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::forward starting" << std::endl;

  // Get vorticity and divergence fields
  auto vorField = fset["air_upward_absolute_vorticity"];
  ASSERT(vorField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(vorField.shape(1) == static_cast<int>(nz_));
  auto divField = fset["air_horizontal_divergence"];
  ASSERT(divField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(divField.shape(1) == static_cast<int>(nz_));

  // Compute stream function and velocity potential
  trans_->inverseLaplacian(vorField);
  trans_->inverseLaplacian(divField);

  // Compute horizontal derivatives
  atlas::Field dPsiDxField;
  atlas::Field dPsiDyField;
  atlas::Field dKhiDxField;
  atlas::Field dKhiDyField;
  trans_->derivative(vorField, dPsiDxField, "x");
  trans_->derivative(vorField, dPsiDyField, "y");
  trans_->derivative(divField, dKhiDxField, "x");
  trans_->derivative(divField, dKhiDyField, "y");

  // Compute u/v
  auto vorView = make_view<double, 2>(vorField);
  auto divView = make_view<double, 2>(divField);
  const auto dPsiDxView = make_view<double, 2>(dPsiDxField);
  const auto dPsiDyView = make_view<double, 2>(dPsiDyField);
  const auto dKhiDxView = make_view<double, 2>(dKhiDxField);
  const auto dKhiDyView = make_view<double, 2>(dKhiDyField);
  for (size_t jz = 0; jz < nz_; ++jz) {
    for (size_t js = 0; js < trans_->ns(); ++js) {
      vorView(js, jz) = dKhiDxView(js, jz) - dPsiDyView(js, jz);
      divView(js, jz) = dKhiDyView(js, jz) + dPsiDxView(js, jz);
    }
  }

  // Reset mean wind profile
  if (jsZero_ >= 0) {
    if (fset.fieldSet().metadata().has("uMeanProfile")
      && fset.fieldSet().metadata().has("vMeanProfile")) {
      const std::vector<double> uMeanProfile =
        fset.fieldSet().metadata().getDoubleVector("uMeanProfile");
      const std::vector<double> vMeanProfile =
        fset.fieldSet().metadata().getDoubleVector("vMeanProfile");
      for (size_t jz = 0; jz < nz_; ++jz) {
        vorView(jsZero_, jz) = uMeanProfile[jz];
        divView(jsZero_, jz) = vMeanProfile[jz];
      }
    }
  }

  // Rename wind fields
  fset["air_upward_absolute_vorticity"].rename("reduced_x_wind");
  fset["air_horizontal_divergence"].rename("reduced_y_wind");

  oops::Log::trace() << classname() << "::forward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::forwardAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::forwardAD starting" << std::endl;

  if (!params_.dipoleTest.value()) {
    // Get vorticity and divergence fields
    auto vorField = fset["reduced_x_wind"];
    ASSERT(vorField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(vorField.shape(1) == static_cast<int>(nz_));
    auto divField = fset["reduced_y_wind"];
    ASSERT(divField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(divField.shape(1) == static_cast<int>(nz_));
    auto vorView = make_view<double, 2>(vorField);
    auto divView = make_view<double, 2>(divField);

    // Save mean wind profile
    if (jsZero_ >= 0) {
      std::vector<double> uMeanProfile(nz_, 0.0);
      std::vector<double> vMeanProfile(nz_, 0.0);
      for (size_t jz = 0; jz < nz_; ++jz) {
        uMeanProfile[jz] = vorView(jsZero_, jz);
        vMeanProfile[jz] = divView(jsZero_, jz);
      }
      fset.fieldSet().metadata().set("uMeanProfile", uMeanProfile);
      fset.fieldSet().metadata().set("vMeanProfile", vMeanProfile);
    }

    // Compute derivatives adjoints
    atlas::Field dPsiDxField;
    atlas::Field dPsiDyField;
    atlas::Field dKhiDxField;
    atlas::Field dKhiDyField;
    trans_->derivative(divField, dPsiDxField, "x", true);
    trans_->derivative(vorField, dPsiDyField, "y", true);
    trans_->derivative(vorField, dKhiDxField, "x", true);
    trans_->derivative(divField, dKhiDyField, "y", true);

    // Compute vorticity and divergence
    const auto dPsiDxView = make_view<double, 2>(dPsiDxField);
    const auto dPsiDyView = make_view<double, 2>(dPsiDyField);
    const auto dKhiDxView = make_view<double, 2>(dKhiDxField);
    const auto dKhiDyView = make_view<double, 2>(dKhiDyField);
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t js = 0; js < trans_->ns(); ++js) {
        vorView(js, jz) = dPsiDxView(js, jz) - dPsiDyView(js, jz);
        divView(js, jz) = dKhiDxView(js, jz) + dKhiDyView(js, jz);
      }
    }

    // Apply inverse Laplacian
    trans_->inverseLaplacian(vorField);
    trans_->inverseLaplacian(divField);
  }

  // Rename wind fields
  fset["reduced_x_wind"].rename("air_upward_absolute_vorticity");
  fset["reduced_y_wind"].rename("air_horizontal_divergence");

  oops::Log::trace() << classname() << "::forwardAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::backward(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::backward starting" << std::endl;

  // Get vorticity and divergence fields
  auto uField = fset["reduced_x_wind"];
  ASSERT(uField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(uField.shape(1) == static_cast<int>(nz_));
  auto vField = fset["reduced_y_wind"];
  ASSERT(vField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(vField.shape(1) == static_cast<int>(nz_));

  // Compute horizontal derivatives
  atlas::Field dUDxField;
  atlas::Field dUDyField;
  atlas::Field dVDxField;
  atlas::Field dVDyField;
  trans_->derivative(uField, dUDxField, "x");
  trans_->derivative(uField, dUDyField, "y");
  trans_->derivative(vField, dVDxField, "x");
  trans_->derivative(vField, dVDyField, "y");

  // Compute vor/div
  auto uView = make_view<double, 2>(uField);
  auto vView = make_view<double, 2>(vField);
  const auto dUDxView = make_view<double, 2>(dUDxField);
  const auto dUDyView = make_view<double, 2>(dUDyField);
  const auto dVDxView = make_view<double, 2>(dVDxField);
  const auto dVDyView = make_view<double, 2>(dVDyField);
  for (size_t jz = 0; jz < nz_; ++jz) {
    for (size_t js = 0; js < trans_->ns(); ++js) {
      uView(js, jz) = dVDxView(js, jz) - dUDyView(js, jz);
      vView(js, jz) = dUDxView(js, jz) + dVDyView(js, jz);
    }
  }

  // Rename wind fields
  fset["reduced_x_wind"].rename("air_upward_absolute_vorticity");
  fset["reduced_y_wind"].rename("air_horizontal_divergence");

  oops::Log::trace() << classname() << "::backward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::backwardAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::backwardAD starting" << std::endl;

  // Get vorticity and divergence fields
  auto uField = fset["air_upward_absolute_vorticity"];
  ASSERT(uField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(uField.shape(1) == static_cast<int>(nz_));
  auto vField = fset["air_horizontal_divergence"];
  ASSERT(vField.shape(0) == static_cast<int>(trans_->ns()));
  ASSERT(vField.shape(1) == static_cast<int>(nz_));

  // Compute horizontal derivatives
  atlas::Field dUDxField;
  atlas::Field dUDyField;
  atlas::Field dVDxField;
  atlas::Field dVDyField;
  trans_->derivative(vField, dUDxField, "x", true);
  trans_->derivative(uField, dUDyField, "y", true);
  trans_->derivative(uField, dVDxField, "x", true);
  trans_->derivative(vField, dVDyField, "y", true);

  // Compute vor/div
  auto uView = make_view<double, 2>(uField);
  auto vView = make_view<double, 2>(vField);
  const auto dUDxView = make_view<double, 2>(dUDxField);
  const auto dUDyView = make_view<double, 2>(dUDyField);
  const auto dVDxView = make_view<double, 2>(dVDxField);
  const auto dVDyView = make_view<double, 2>(dVDyField);
  for (size_t jz = 0; jz < nz_; ++jz) {
    for (size_t js = 0; js < trans_->ns(); ++js) {
      uView(js, jz) = dUDxView(js, jz) - dUDyView(js, jz);
      vView(js, jz) = dVDxView(js, jz) + dVDyView(js, jz);
    }
  }

  // Rename wind fields
  fset["air_upward_absolute_vorticity"].rename("reduced_x_wind");
  fset["air_horizontal_divergence"].rename("reduced_y_wind");

  oops::Log::trace() << classname() << "::backwardAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorDivToRedWind::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
