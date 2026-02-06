/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierSpectralVorDivToGridWind.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierSpectralVorDivToGridWind>
  makerBifourierSpectralVorDivToGridWind_("BifourierSpectralVorDivToGridWind");

// -----------------------------------------------------------------------------

BifourierSpectralVorDivToGridWind::BifourierSpectralVorDivToGridWind(
  const oops::GeometryData & outerGeometryData,
  const oops::Variables & outerVars,
  const eckit::Configuration &,
  const Parameters_ & params,
  const oops::FieldSet3D & xb,
  const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    comm_(outerGeometryData.comm()),
    outerVars_(outerVars),
    params_(params),
    fftBackend_(params_.transform.value().fftBackend.value()),
    data_(xb.validTime(), comm_)
{
  oops::Log::trace() << classname() << "::BifourierSpectralVorDivToGridWind starting"
    << std::endl;

  // Inner variables
  if (!params_.backward.value()) {
    // Wind variables
    oops::Variables windVars;
    if (params_.outerSphericalWinds.value()) {
      windVars.push_back(outerVars["eastward_wind"]);
      windVars.push_back(outerVars["northward_wind"]);
    } else {
      windVars.push_back(outerVars["geographical_x_wind"]);
      windVars.push_back(outerVars["geographical_y_wind"]);
    }

    // Add vor/div to inner variables
    innerVars_.push_back("air_upward_absolute_vorticity");
    innerVars_.push_back("air_horizontal_divergence");

    // Number of levels
    nz_ = windVars[0].getLevels();
    innerVars_["air_upward_absolute_vorticity"].setLevels(nz_);
    innerVars_["air_horizontal_divergence"].setLevels(nz_);

    // Add other variables
    for (const auto & var : outerVars) {
      if (!windVars.has(var)) {
        innerVars_.push_back(var);
      }
    }
  } else {
    // Wind variables
    oops::Variables windVars;
    if (params_.outerSphericalWinds.value()) {
      windVars.push_back("eastward_wind");
      windVars.push_back("northward_wind");
    } else {
      windVars.push_back("geographical_x_wind");
      windVars.push_back("geographical_y_wind");
    }

    // Number of levels
    nz_ = outerVars["air_upward_absolute_vorticity"].getLevels();
    windVars[0].setLevels(nz_);
    windVars[1].setLevels(nz_);

    // Add wind to inner variables
    innerVars_.push_back(windVars[0]);
    innerVars_.push_back(windVars[1]);

    // Add other variables
    for (const auto & var : outerVars) {
      if (!((var.name() == "air_upward_absolute_vorticity")
        || (var.name() == "air_horizontal_divergence"))) {
        innerVars_.push_back(var);
      }
    }
  }

  if (!params_.backward.value()) {
    // Create spectral transform
    trans_ = transStore_.setupTransform(outerGeometryData, innerVars_,
      params_.transform.value());

    // Create inner GeometryData
    innerGeometryData_.reset(new oops::GeometryData(trans_->spFspace(),
      outerGeometryData.fieldSet(), outerGeometryData.levelsAreTopDown(),
      outerGeometryData.comm()));
  } else {
    // Retrieve spectral transform
    trans_ = transStore_.retrieveTransform(outerGeometryData);

    // Set inner GeometryData
    // TODO(Benjamin): avoid this
    innerGeometryData_.reset(new oops::GeometryData(trans_->geometryData().functionSpace(),
      trans_->geometryData().fieldSet(), trans_->geometryData().levelsAreTopDown(),
      trans_->geometryData().comm()));
  }

  // Prepare biperiodization if needed
  const auto &biperParams = params_.biperParams.value();
  oops::Variables biperVars;

  // Get grid function space
  const atlas::functionspace::StructuredColumns fs(trans_->geometryData().functionSpace());

  // Get lon/lat view
  const auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());

  if (fftBackend_ == "fftw") {
    // Get js for (jk, jl) = (0, 0)
    jsZero_ = -1;
    for (size_t js = 0; js < trans_->ns(); ++js) {
      if ((trans_->jk(js) == 0) && (trans_->jl(js) == 0)) {
        jsZero_ = js;
      }
    }

    if (biperParams != boost::none) {
      // Check biperiodization parameters (should not change the grid size)
      ASSERT(biperParams->innerExtNx.value() == biperParams->outerExtNx.value());
      ASSERT(biperParams->innerExtNy.value() == biperParams->outerExtNy.value());
    }

    // Create map factor field
    atlas::Field mapFactorField = fs.createField<double>(
      atlas::option::name("map_factor") | atlas::option::levels(1));

    // Add map factor to FieldSet
    data_.add(mapFactorField);

    // Get map factor view
    auto mapFactorView = make_view<double, 2>(mapFactorField);

    // Degree to radian
    const double deg2rad = M_PI/180.0;

    // Get projection parameter
    const eckit::LocalConfiguration projConf = fs.grid().projection().spec();
    const double latitude0 = projConf.getDouble("latitude0")*deg2rad;

    // Get pole
    double pole;
    if (latitude0 == 0.0) {
      pole = 0.0;
    } else if (latitude0 > 0.0) {
      pole = 1.0;
    } else {
      pole = -1.0;
    }

    // Compute map factor constants
    const double sina = pole*std::sin(latitude0);
    const double rho0 = std::abs(sina) > 0.0 ?
       std::pow(std::cos(latitude0), 1.0-sina)*std::pow(1.0+sina, sina)/sina : 1.0;

    // Compute map factor
    for (int jnode = 0; jnode < mapFactorField.shape(0); ++jnode) {
      const double latRad = lonlatView(jnode, 1)*deg2rad;
      mapFactorView(jnode, 0) = rho0*sina/std::cos(latRad)
        *std::pow(std::tan((0.25*M_PI)-(pole*0.5*latRad)), sina);
    }

    if (biperParams != boost::none) {
      // Add biperiodization variable
      biperVars.push_back("map_factor");
      biperVars["map_factor"].setLevels(1);
    }
  }

  if (params_.outerSphericalWinds.value()) {
    // Get eastward and northward winds from grid winds

    // Create coefficients fields
    atlas::Field dxDlonField = fs.createField<double>(
      atlas::option::name("dxDlon") | atlas::option::levels(1));
    atlas::Field dxDlatField = fs.createField<double>(
      atlas::option::name("dxDlat") | atlas::option::levels(1));
    atlas::Field dyDlonField = fs.createField<double>(
      atlas::option::name("dyDlon") | atlas::option::levels(1));
    atlas::Field dyDlatField = fs.createField<double>(
      atlas::option::name("dyDlat") | atlas::option::levels(1));

    // Add coefficients to FieldSet
    data_.add(dxDlonField);
    data_.add(dxDlatField);
    data_.add(dyDlonField);
    data_.add(dyDlatField);

    // Get coefficients views
    auto dxDlonView = make_view<double, 2>(dxDlonField);
    auto dxDlatView = make_view<double, 2>(dxDlatField);
    auto dyDlonView = make_view<double, 2>(dyDlonField);
    auto dyDlatView = make_view<double, 2>(dyDlatField);

    for (int jnode = 0; jnode < dxDlonField.shape(0); ++jnode) {
      // Get local point
      atlas::PointLonLat p({lonlatView(jnode, 0), lonlatView(jnode, 1)});

      // Get local Jacobian
      dxDlonView(jnode, 0) = fs.grid().projection().jacobian(p).dx_dlon();
      dxDlatView(jnode, 0) = fs.grid().projection().jacobian(p).dx_dlat();
      dyDlonView(jnode, 0) = fs.grid().projection().jacobian(p).dy_dlon();
      dyDlatView(jnode, 0) = fs.grid().projection().jacobian(p).dy_dlat();

      // Normalize Jacobian
      const double dlonNorm = 1.0/std::sqrt(
        dxDlonView(jnode, 0)*dxDlonView(jnode, 0)+dyDlonView(jnode, 0)*dyDlonView(jnode, 0));
      const double dlatNorm = 1.0/std::sqrt(
        dxDlatView(jnode, 0)*dxDlatView(jnode, 0)+dyDlatView(jnode, 0)*dyDlatView(jnode, 0));
      dxDlonView(jnode, 0) *= dlonNorm;
      dxDlatView(jnode, 0) *= dlatNorm;
      dyDlonView(jnode, 0) *= dlonNorm;
      dyDlatView(jnode, 0) *= dlatNorm;
    }

    if (biperParams != boost::none) {
      // Add biperiodization variables
      biperVars.push_back("dxDlon");
      biperVars["dxDlon"].setLevels(1);
      biperVars.push_back("dxDlat");
      biperVars["dxDlat"].setLevels(1);
      biperVars.push_back("dyDlon");
      biperVars["dyDlon"].setLevels(1);
      biperVars.push_back("dyDlat");
      biperVars["dyDlat"].setLevels(1);
    }
  }

  if (biperParams != boost::none) {
    // Setup biperiodization implementation
    BiperiodizationImpl biper(outerGeometryData, biperVars, *biperParams);

    // Apply biperiodization leftInverseMultiply to go to biperiodization inner geometry
    biper.leftInverseMultiply(data_.fieldSet());

    // Apply biperiodization multiply
    biper.multiply(data_.fieldSet());
  }

  oops::Log::trace() << classname() << "::BifourierSpectralVorDivToGridWind done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Forward application
    forward(fset, innerVars_);
  } else {
    // Backward application
    backward(fset, outerVars_);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (!params_.backward.value()) {
    // Forward application, adjoint
    forwardAD(fset, innerVars_);
  } else {
    // Backward application, adjoint
    backwardAD(fset, outerVars_);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Backward application
    backward(fset, innerVars_);
  } else {
    // Forward application
    forward(fset, outerVars_);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::forward(oops::FieldSet3D & fset,
                                                const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::forward starting" << std::endl;

  if (fftBackend_ == "fftw") {
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
  } else if (fftBackend_ == "ectrans") {
    // Set number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", nz_);
  }

  // Create temporary fieldset with the same metadata
  atlas::FieldSet fsetTmp;
  fsetTmp.metadata() = fset.fieldSet().metadata();

  // Copy FieldSet
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp, vars);

  // Remove variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), vars.variables());

  // Inverse spectral transform
  trans_->sp2gp(fsetTmp, fset.fieldSet(), vars);

  // Rename wind fields
  fset["air_upward_absolute_vorticity"].rename("reduced_x_wind");
  fset["air_horizontal_divergence"].rename("reduced_y_wind");

  // Get reduced wind fields
  auto uField = fset["reduced_x_wind"];
  auto vField = fset["reduced_y_wind"];

  // Get views
  auto uView = atlas::array::make_view<double, 2>(uField);
  auto vView = atlas::array::make_view<double, 2>(vField);

  if (fftBackend_ == "fftw") {
    // Get map factor view
    const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

    // Multiply by the map factor
    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        uView(jnode, jlevel) *= mapFactorView(jnode, 0);
        vView(jnode, jlevel) *= mapFactorView(jnode, 0);
      }
    }
  }

  if (params_.outerSphericalWinds.value()) {
    // Get coefficients views
    const auto dxDlonView = atlas::array::make_view<double, 2>(data_["dxDlon"]);
    const auto dxDlatView = atlas::array::make_view<double, 2>(data_["dxDlat"]);
    const auto dyDlonView = atlas::array::make_view<double, 2>(data_["dyDlon"]);
    const auto dyDlatView = atlas::array::make_view<double, 2>(data_["dyDlat"]);

    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        const double uSave = uView(jnode, jlevel);
        const double vSave = vView(jnode, jlevel);
        uView(jnode, jlevel) = uSave*dxDlonView(jnode, 0) + vSave*dyDlonView(jnode, 0);
        vView(jnode, jlevel) = uSave*dxDlatView(jnode, 0) + vSave*dyDlatView(jnode, 0);
      }
    }

    // Rename wind fields
    fset["reduced_x_wind"].rename("eastward_wind");
    fset["reduced_y_wind"].rename("northward_wind");
  } else {
    // Rename wind fields
    fset["reduced_x_wind"].rename("geographical_x_wind");
    fset["reduced_y_wind"].rename("geographical_y_wind");
  }

  if (fftBackend_ == "ectrans") {
    // Reset number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", 0);
  }

  oops::Log::trace() << classname() << "::forward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::forwardAD(oops::FieldSet3D & fset,
                                                  const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::forwardAD starting" << std::endl;

  if (params_.outerSphericalWinds.value()) {
    // Rename wind fields
    fset["eastward_wind"].rename("reduced_x_wind");
    fset["northward_wind"].rename("reduced_y_wind");
  } else {
    // Rename wind fields
    fset["geographical_x_wind"].rename("reduced_x_wind");
    fset["geographical_y_wind"].rename("reduced_y_wind");
  }

  // Get reduced wind fields
  auto uField = fset["reduced_x_wind"];
  auto vField = fset["reduced_y_wind"];

  // Get views
  auto uView = atlas::array::make_view<double, 2>(uField);
  auto vView = atlas::array::make_view<double, 2>(vField);

  if (params_.outerSphericalWinds.value() && !params_.dipoleTest.value()) {
    // Get coefficients views
    const auto dxDlonView = atlas::array::make_view<double, 2>(data_["dxDlon"]);
    const auto dxDlatView = atlas::array::make_view<double, 2>(data_["dxDlat"]);
    const auto dyDlonView = atlas::array::make_view<double, 2>(data_["dyDlon"]);
    const auto dyDlatView = atlas::array::make_view<double, 2>(data_["dyDlat"]);

    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        const double uSave = uView(jnode, jlevel);
        const double vSave = vView(jnode, jlevel);
        uView(jnode, jlevel) = uSave*dxDlonView(jnode, 0) + vSave*dxDlatView(jnode, 0);
        vView(jnode, jlevel) = uSave*dyDlonView(jnode, 0) + vSave*dyDlatView(jnode, 0);
      }
    }
  }

  if (fftBackend_ == "fftw" && !params_.dipoleTest.value()) {
    // Get map factor view
    const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

    // Multiply by the map factor
    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        uView(jnode, jlevel) *= mapFactorView(jnode, 0);
        vView(jnode, jlevel) *= mapFactorView(jnode, 0);
      }
    }
  }

  // Rename wind fields
  fset["reduced_x_wind"].rename("air_upward_absolute_vorticity");
  fset["reduced_y_wind"].rename("air_horizontal_divergence");

  if (fftBackend_ == "ectrans" && !params_.dipoleTest.value()) {
    // Set number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", nz_);
  }

  // Create temporary fieldset with the same metadata
  atlas::FieldSet fsetTmp;
  fsetTmp.metadata() = fset.fieldSet().metadata();

  // Direct spectral transform
  trans_->sp2gpAdj(fset.fieldSet(), fsetTmp, vars);

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), vars.variables());

  // Copy FieldSet
  trans_->copyFieldSet(fsetTmp, fset.fieldSet(), vars);

  if (fftBackend_ == "fftw" && !params_.dipoleTest.value()) {
    // Get vorticity and divergence fields
    auto vorField = fset["air_upward_absolute_vorticity"];
    ASSERT(vorField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(vorField.shape(1) == static_cast<int>(nz_));
    auto divField = fset["air_horizontal_divergence"];
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
  } else if (fftBackend_ == "ectrans" && !params_.dipoleTest.value()) {
    // Reset number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", 0);
  }

  oops::Log::trace() << classname() << "::forwardAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::backward(oops::FieldSet3D & fset,
                                                 const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::backward starting" << std::endl;

  if (params_.outerSphericalWinds.value()) {
    // Rename wind fields
    fset["eastward_wind"].rename("reduced_x_wind");
    fset["northward_wind"].rename("reduced_y_wind");
  } else {
    // Rename wind fields
    fset["geographical_x_wind"].rename("reduced_x_wind");
    fset["geographical_y_wind"].rename("reduced_y_wind");
  }

  // Get reduced wind fields
  auto uField = fset["reduced_x_wind"];
  auto vField = fset["reduced_y_wind"];

  // Get views
  auto uView = atlas::array::make_view<double, 2>(uField);
  auto vView = atlas::array::make_view<double, 2>(vField);

  if (params_.outerSphericalWinds.value()) {
    // Get coefficients views
    const auto dxDlonView = atlas::array::make_view<double, 2>(data_["dxDlon"]);
    const auto dxDlatView = atlas::array::make_view<double, 2>(data_["dxDlat"]);
    const auto dyDlonView = atlas::array::make_view<double, 2>(data_["dyDlon"]);
    const auto dyDlatView = atlas::array::make_view<double, 2>(data_["dyDlat"]);

    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      // Inverse matrix
      const double det = dxDlonView(jnode, 0)*dyDlatView(jnode, 0)
        -dyDlonView(jnode, 0)*dxDlatView(jnode, 0);
      const double dlonDx = dyDlatView(jnode, 0)/det;
      const double dlonDy = -dyDlonView(jnode, 0)/det;
      const double dlatDx = -dxDlatView(jnode, 0)/det;
      const double dlatDy = dxDlonView(jnode, 0)/det;

      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        // Apply inverse matrix
        const double uSave = uView(jnode, jlevel);
        const double vSave = vView(jnode, jlevel);
        uView(jnode, jlevel) = uSave*dlonDx + vSave*dlonDy;
        vView(jnode, jlevel) = uSave*dlatDx + vSave*dlatDy;
      }
    }
  }

  if (fftBackend_ == "fftw") {
    // Get map factor view
    const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

    // Divide by the map factor
    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        uView(jnode, jlevel) /= mapFactorView(jnode, 0);
        vView(jnode, jlevel) /= mapFactorView(jnode, 0);
      }
    }
  }

  // Rename wind fields
  fset["reduced_x_wind"].rename("air_upward_absolute_vorticity");
  fset["reduced_y_wind"].rename("air_horizontal_divergence");

  if (fftBackend_ == "ectrans") {
    // Set number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", nz_);
  }

  // Create temporary fieldset with the same metadata
  atlas::FieldSet fsetTmp;
  fsetTmp.metadata() = fset.fieldSet().metadata();

  // Direct spectral transform
  trans_->gp2sp(fset.fieldSet(), fsetTmp, vars);

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), vars.variables());

  // Copy FieldSet
  trans_->copyFieldSet(fsetTmp, fset.fieldSet(), vars);

  if (fftBackend_ == "fftw") {
    // Get vorticity and divergence fields
    auto vorField = fset["air_upward_absolute_vorticity"];
    ASSERT(vorField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(vorField.shape(1) == static_cast<int>(nz_));
    auto divField = fset["air_horizontal_divergence"];
    ASSERT(divField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(divField.shape(1) == static_cast<int>(nz_));

    // Compute horizontal derivatives
    atlas::Field dUDxField;
    atlas::Field dUDyField;
    atlas::Field dVDxField;
    atlas::Field dVDyField;
    trans_->derivative(vorField, dUDxField, "x");
    trans_->derivative(vorField, dUDyField, "y");
    trans_->derivative(divField, dVDxField, "x");
    trans_->derivative(divField, dVDyField, "y");

    // Compute vor/div
    auto vorView = make_view<double, 2>(vorField);
    auto divView = make_view<double, 2>(divField);
    const auto dUDxView = make_view<double, 2>(dUDxField);
    const auto dUDyView = make_view<double, 2>(dUDyField);
    const auto dVDxView = make_view<double, 2>(dVDxField);
    const auto dVDyView = make_view<double, 2>(dVDyField);
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t js = 0; js < trans_->ns(); ++js) {
        vorView(js, jz) = dVDxView(js, jz) - dUDyView(js, jz);
        divView(js, jz) = dUDxView(js, jz) + dVDyView(js, jz);
      }
    }
  } else if (fftBackend_ == "ectrans") {
    // Reset number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", 0);
  }

  oops::Log::trace() << classname() << "::backward done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::backwardAD(oops::FieldSet3D & fset,
                                                   const oops::Variables & vars) const {
  oops::Log::trace() << classname() << "::backwardAD starting" << std::endl;

  if (fftBackend_ == "fftw") {
    // Get vorticity and divergence fields
    auto vorField = fset["air_upward_absolute_vorticity"];
    ASSERT(vorField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(vorField.shape(1) == static_cast<int>(nz_));
    auto divField = fset["air_horizontal_divergence"];
    ASSERT(divField.shape(0) == static_cast<int>(trans_->ns()));
    ASSERT(divField.shape(1) == static_cast<int>(nz_));

    // Compute horizontal derivatives
    atlas::Field dUDxField;
    atlas::Field dUDyField;
    atlas::Field dVDxField;
    atlas::Field dVDyField;
    trans_->derivative(divField, dUDxField, "x", true);
    trans_->derivative(vorField, dUDyField, "y", true);
    trans_->derivative(vorField, dVDxField, "x", true);
    trans_->derivative(divField, dVDyField, "y", true);

    // Compute vor/div
    auto vorView = make_view<double, 2>(vorField);
    auto divView = make_view<double, 2>(divField);
    const auto dUDxView = make_view<double, 2>(dUDxField);
    const auto dUDyView = make_view<double, 2>(dUDyField);
    const auto dVDxView = make_view<double, 2>(dVDxField);
    const auto dVDyView = make_view<double, 2>(dVDyField);
    for (size_t jz = 0; jz < nz_; ++jz) {
      for (size_t js = 0; js < trans_->ns(); ++js) {
        vorView(js, jz) = dUDxView(js, jz) - dUDyView(js, jz);
        divView(js, jz) = dVDxView(js, jz) + dVDyView(js, jz);
      }
    }
  } else if (fftBackend_ == "ectrans") {
    // Set number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", nz_);
  }

  // Create temporary fieldset with the same metadata
  atlas::FieldSet fsetTmp;
  fsetTmp.metadata() = fset.fieldSet().metadata();

  // Copy FieldSet
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp, vars);

  // Remove variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), vars.variables());

  // Inverse spectral transform
  trans_->gp2spAdj(fsetTmp, fset.fieldSet(), vars);

  // Rename wind fields
  fset["air_upward_absolute_vorticity"].rename("reduced_x_wind");
  fset["air_horizontal_divergence"].rename("reduced_y_wind");

  // Get reduced wind fields
  auto uField = fset["reduced_x_wind"];
  auto vField = fset["reduced_y_wind"];

  // Get views
  auto uView = atlas::array::make_view<double, 2>(uField);
  auto vView = atlas::array::make_view<double, 2>(vField);

  if (fftBackend_ == "fftw") {
    // Get map factor view
    const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

    // Divide by the map factor
    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        uView(jnode, jlevel) /= mapFactorView(jnode, 0);
        vView(jnode, jlevel) /= mapFactorView(jnode, 0);
      }
    }
  }

  if (params_.outerSphericalWinds.value()) {
    // Get coefficients views
    const auto dxDlonView = atlas::array::make_view<double, 2>(data_["dxDlon"]);
    const auto dxDlatView = atlas::array::make_view<double, 2>(data_["dxDlat"]);
    const auto dyDlonView = atlas::array::make_view<double, 2>(data_["dyDlon"]);
    const auto dyDlatView = atlas::array::make_view<double, 2>(data_["dyDlat"]);

    for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
      // Inverse matrix
      const double det = dxDlonView(jnode, 0)*dyDlatView(jnode, 0)
        -dyDlonView(jnode, 0)*dxDlatView(jnode, 0);
      const double dlonDx = dyDlatView(jnode, 0)/det;
      const double dlonDy = -dyDlonView(jnode, 0)/det;
      const double dlatDx = -dxDlatView(jnode, 0)/det;
      const double dlatDy = dxDlonView(jnode, 0)/det;

      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        // Apply inverse matrix
        const double uSave = uView(jnode, jlevel);
        const double vSave = vView(jnode, jlevel);
        uView(jnode, jlevel) = uSave*dlonDx + vSave*dlatDx;
        vView(jnode, jlevel) = uSave*dlonDy + vSave*dlatDy;
      }
    }
  }

  if (params_.outerSphericalWinds.value()) {
    // Rename wind fields
    fset["reduced_x_wind"].rename("eastward_wind");
    fset["reduced_y_wind"].rename("northward_wind");
  } else {
    // Rename wind fields
    fset["reduced_x_wind"].rename("geographical_x_wind");
    fset["reduced_y_wind"].rename("geographical_y_wind");
  }

  if (fftBackend_ == "ectrans") {
    // Reset number of levels for the wind transform
    fset.fieldSet().metadata().set("nvordiv", 0);
  }

  oops::Log::trace() << classname() << "::backwardAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralVorDivToGridWind::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
