/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/RedWindToGeoWind.h"

#include "atlas/field.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<RedWindToGeoWind>
  makerRedWindToGeoWind_("RedWindToGeoWind");

// -----------------------------------------------------------------------------

RedWindToGeoWind::RedWindToGeoWind(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConfig,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params),
    data_(xb.validTime(), comm_)
{
  oops::Log::trace() << classname() << "::RedWindToGeoWind starting" << std::endl;

  // Process variables
  size_t nz;
  if (params_.outerSphericalWinds.value()) {
    // Check variables
    ASSERT(innerVars_.has("eastward_wind") && innerVars_.has("northward_wind"));

    // Get number of levels
    nz = innerVars_["eastward_wind"].getLevels();

    // Remove spherical winds from inner variables
    innerVars_ -= innerVars_["eastward_wind"];
    innerVars_ -= innerVars_["northward_wind"];
  } else {
    // Check variables
    ASSERT(innerVars_.has("geographical_x_wind") && innerVars_.has("geographical_y_wind"));

    // Get number of levels
    nz = innerVars_["geographical_x_wind"].getLevels();

    // Remove grid winds from inner variables
    innerVars_ -= innerVars_["geographical_x_wind"];
    innerVars_ -= innerVars_["geographical_y_wind"];
  }

  // Add reduced wind to inner variables
  innerVars_.push_back("reduced_x_wind");
  innerVars_["reduced_x_wind"].setLevels(nz);
  innerVars_.push_back("reduced_y_wind");
  innerVars_["reduced_y_wind"].setLevels(nz);

  // Prepare biperiodization if needed
  const auto &biperParams = params_.biperParams.value();
  if (biperParams != boost::none) {
    // Check biperiodization parameters (should not change the grid size)
    ASSERT(biperParams->innerExtNx.value() == biperParams->outerExtNx.value());
    ASSERT(biperParams->innerExtNy.value() == biperParams->outerExtNy.value());
  }
  oops::Variables biperVars;

  // Get grid function space
  const atlas::functionspace::StructuredColumns fs(innerGeometryData_.functionSpace());

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

  // Get lon/lat view
  const auto lonlatView = atlas::array::make_view<double, 2>(fs.lonlat());

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

  oops::Log::trace() << classname() << "::RedWindToGeoWind done" << std::endl;
}

// -----------------------------------------------------------------------------

void RedWindToGeoWind::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Get map factor view
  const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

  // Get reduced wind fields
  auto uField = fset["reduced_x_wind"];
  auto vField = fset["reduced_y_wind"];

  // Get views
  auto uView = atlas::array::make_view<double, 2>(uField);
  auto vView = atlas::array::make_view<double, 2>(vField);

  // Multiply by the map factor
  for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
    for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
      uView(jnode, jlevel) *= mapFactorView(jnode, 0);
      vView(jnode, jlevel) *= mapFactorView(jnode, 0);
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

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void RedWindToGeoWind::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

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
      for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
        const double uSave = uView(jnode, jlevel);
        const double vSave = vView(jnode, jlevel);
        uView(jnode, jlevel) = uSave*dxDlonView(jnode, 0) + vSave*dxDlatView(jnode, 0);
        vView(jnode, jlevel) = uSave*dyDlonView(jnode, 0) + vSave*dyDlatView(jnode, 0);
      }
    }
  }

  // Get map factor view
  const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

  // Multiply by the map factor
  for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
    for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
      uView(jnode, jlevel) *= mapFactorView(jnode, 0);
      vView(jnode, jlevel) *= mapFactorView(jnode, 0);
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void RedWindToGeoWind::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

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

  // Get map factor view
  const auto mapFactorView = atlas::array::make_view<double, 2>(data_["map_factor"]);

  // Multiply by the map factor
  for (int jnode = 0; jnode < uField.shape(0); ++jnode) {
    for (int jlevel = 0; jlevel < uField.shape(1); ++jlevel) {
      uView(jnode, jlevel) /= mapFactorView(jnode, 0);
      vView(jnode, jlevel) /= mapFactorView(jnode, 0);
    }
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>>
  RedWindToGeoWind::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> pairs;
  if (params_.outputFile.value() != boost::none) {
    pairs.push_back(std::pair<eckit::LocalConfiguration, oops::FieldSet3D>(
      *params_.outputFile.value(), data_));
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}
// -----------------------------------------------------------------------------

void RedWindToGeoWind::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
