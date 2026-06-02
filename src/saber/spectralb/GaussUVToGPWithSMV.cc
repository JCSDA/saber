/*
 * (C) Crown Copyright 2022-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "atlas/util/Constants.h"

#include "eckit/exception/Exceptions.h"

#include "mo/eval_cloud_ice_mixing_ratio.h"
#include "mo/eval_cloud_liquid_water_mixing_ratio.h"
#include "mo/eval_dry_air_density.h"
#include "mo/eval_total_water_mixing_ratio.h"
#include "mo/eval_water_vapor_mixing_ratio.h"

#include "oops/util/for_each.h"

#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/interpolation/SMVInterpWrapper.h"
#include "saber/oops/Utilities.h"
#include "saber/spectralb/GaussUVToGPWithSMV.h"


#define ERR(e)                                      \
  {                                                 \
    throw eckit::Exception(nc_strerror(e), Here()); \
  }

namespace saber {
namespace spectralb {

namespace {

double normfield(const eckit::mpi::Comm& comm, const atlas::Field& fld) {
  auto view = atlas::array::make_view<const double, 2>(fld);

  std::size_t i(0);
  double zz(0.0);

  const atlas::idx_t sizeOwned = fld.functionspace().type() == "PointCloud"
                                     ? fld.functionspace().size()
                                     : util::getSizeOwned(fld.functionspace());

  for (atlas::idx_t jn = 0; jn < sizeOwned; ++jn) {
    for (atlas::idx_t jl = 0; jl < fld.shape(1); ++jl, ++i) {
      zz += view(jn, jl) * view(jn, jl);
    }
  }
  comm.allReduceInPlace(i, eckit::mpi::sum());
  comm.allReduceInPlace(zz, eckit::mpi::sum());

  return std::sqrt(zz) / static_cast<double>(i);
}

// -----------------------------------------------------------------------------

atlas::Field allocateRHSVec(const atlas::FunctionSpace& gaussFS,
                            const std::size_t& levels) {
  atlas::Field rhsvec = gaussFS.createField<double>(
      atlas::option::name("rhsvec") | atlas::option::variables(2) |
      atlas::option::levels(levels));
  auto rhsVecView = atlas::array::make_view<double, 3>(rhsvec);
  rhsVecView.assign(0.0);
  return rhsvec;
}

// -----------------------------------------------------------------------------

void populateRHSVec(const atlas::Field& rhoState, const atlas::Field& coriolis,
                    const atlas::FieldSet& fset, atlas::Field& rhsvec) {
  auto uView = atlas::array::make_view<const double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<const double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto coriolisView = atlas::array::make_view<const double, 2>(coriolis);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  util::for_each_index(
      util::make_index_space_2d(util::IndexRange::exclude_halo, rhoState),
      [=](const atlas::idx_t jn, const atlas::idx_t jl) mutable {
        rhsvecView(jn, jl, atlas::LON) =
            -coriolisView(jn, 0) * rhoStateView(jn, jl) * vView(jn, jl);
        rhsvecView(jn, jl, atlas::LAT) =
            coriolisView(jn, 0) * rhoStateView(jn, jl) * uView(jn, jl);
      });
  rhsvec.set_dirty();
}

// -----------------------------------------------------------------------------

void populateRHSVecAdj(const atlas::Field& rhoState,
                       const atlas::Field& coriolis, atlas::FieldSet& fset,
                       atlas::Field& rhsvec) {
  auto uView = atlas::array::make_view<double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto coriolisView = atlas::array::make_view<const double, 2>(coriolis);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  util::for_each_index(
      util::make_index_space_2d(util::IndexRange::exclude_halo, rhoState),
      [=](const atlas::idx_t jn, const atlas::idx_t jl) mutable {
        vView(jn, jl) += -coriolisView(jn, 0) * rhoStateView(jn, jl) *
                         rhsvecView(jn, jl, atlas::LON);
        uView(jn, jl) += coriolisView(jn, 0) * rhoStateView(jn, jl) *
                         rhsvecView(jn, jl, atlas::LAT);
        rhsvecView(jn, jl, atlas::LON) = 0.0;
        rhsvecView(jn, jl, atlas::LAT) = 0.0;
      });
  rhsvec.set_dirty();
  fset["eastward_wind"].set_dirty();
  fset["northward_wind"].set_dirty();
}

// -----------------------------------------------------------------------------

atlas::FieldSet populateFields(const atlas::FieldSet& geomfields,
                               const atlas::FieldSet& statefields) {
  // Need to setup derived state fields that we need (done on model grid).
  const std::vector<std::string> requiredStateVariables{
      "air_pressure_levels_minus_one",
      "dry_air_density_levels_minus_one",
      "height_above_mean_sea_level",
      "height_above_mean_sea_level_levels",
      "cloud_ice_mixing_ratio_wrt_dry_air",
      "cloud_liquid_water_mixing_ratio_wrt_dry_air",
      "rain_mixing_ratio_wrt_dry_air",
      "water_vapor_mixing_ratio_wrt_dry_air",
      "total_water_mixing_ratio_wrt_dry_air",
      "air_potential_temperature",
      "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water",
      "cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water",
      "cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water"};

  const std::vector<std::string> requiredGeometryVariables{
      "height_above_mean_sea_level_levels", "height_above_mean_sea_level"};

  const std::string outputVariable{"dry_air_density_levels_minus_one"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  for (auto& s : requiredStateVariables) {
    if (!statefields.has(s)) {
      oops::Log::error() << "::gaussuvtogp::populateFields variable " << s
                         << " is not part of state object." << std::endl;
    }
  }

  atlas::FieldSet tempfields;

  for (const auto& s : requiredStateVariables) {
    tempfields.add(statefields[s]);
  }

  for (const auto& s : requiredGeometryVariables) {
    tempfields.add(geomfields[s]);
  }

  mo::eval_total_water_mixing_ratio_wrt_dry_air_nl(tempfields);
  mo::eval_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
      tempfields);
  mo::eval_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
      tempfields);
  mo::eval_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_nl(
      tempfields);
  mo::eval_dry_air_density_from_pressure_levels_minus_one_nl(tempfields);

  return atlas::FieldSet(tempfields[outputVariable]);
}

// -----------------------------------------------------------------------------

void interpolateCSToGauss(const oops::GeometryData& outerGeometryData,
                          const atlas::FieldSet& CSFields,
                          atlas::FieldSet& gaussFields,
                          const bool normaliseInterpolation) {
  const std::string startComm = eckit::mpi::comm().name();
  eckit::mpi::setCommDefault(outerGeometryData.comm().name());

  // Set up function spaces for interpolation
  std::string s("dry_air_density_levels_minus_one");
  const auto CSFunctionSpace =
      atlas::functionspace::NodeColumns(CSFields[s].functionspace());
  const auto gaussFunctionSpace = atlas::functionspace::StructuredColumns(
      outerGeometryData.functionSpace());
  const atlas::StructuredGrid gaussGrid = gaussFunctionSpace.grid();

  // Do interpolation
  interpolation::SMVInterpWrapper interpolator =
      interpolation::SMVInterpWrapper(CSFunctionSpace, gaussFunctionSpace,
                                      gaussGrid, false, normaliseInterpolation);

  const std::vector<atlas::util::Config> fieldConfigs = {
      atlas::option::name(s) | atlas::option::levels(CSFields[s].shape(1)) |
      atlas::option::halo(gaussFunctionSpace.halo())};
  interpolator.doInterpolation(fieldConfigs, CSFields, gaussFields);

  eckit::mpi::setCommDefault(startComm);
}

// -----------------------------------------------------------------------------

atlas::Field createCoriolis(const atlas::Field& scstate) {
  // create coriolis parameter.
  // need to get the latitude at each grid point from the mesh.
  static constexpr double twoOmega = 2.0 * 7.292116E-5;
  atlas::Field coriolis = scstate.functionspace().createField<double>(
      atlas::option::name("coriolis") | atlas::option::levels(1));

  auto coriolisView = atlas::array::make_view<double, 2>(coriolis);
  auto sc = atlas::functionspace::StructuredColumns(scstate.functionspace());

  atlas_omp_parallel_for(atlas::idx_t j = sc.j_begin(); j < sc.j_end(); ++j) {
    for (atlas::idx_t i = sc.i_begin(j); i < sc.i_end(j); ++i) {
      auto jn = sc.index(i, j);
      coriolisView(jn, 0) =
          twoOmega * sin(atlas::util::Constants::degreesToRadians() *
                         sc.grid().lonlat(i, j).lat());
    }
  }

  return coriolis;
}

// -----------------------------------------------------------------------------

atlas::FieldSet createAugmentedState(
    const oops::GeometryData& outerGeometryData,
    const atlas::FieldSet& xb,
    const bool normaliseInterpolation) {
  const std::string error_message(
      "ERROR - saber block gauss winds to geostrophic pressure with SMV failed.");
  const std::string s("dry_air_density_levels_minus_one");
  const double tolerance(1.0e-6);

  // Error trap that xb has field "s" even if it is not populated.
  if (!xb.has(s)) {
    oops::Log::error() << "Expect " << s << " to be allocated in xb fieldset."
                       << std::endl;
    throw eckit::Exception(error_message, Here());
  }

  // Ensure comm is set to the outerGeometryData comm.
  eckit::mpi::setCommDefault(outerGeometryData.comm().name());

  // Get and check model grid name
  enum ModelGrid { gauss, cubedsphere, unknown };
  ModelGrid modelGrid = unknown;
  if (xb[s].functionspace().type() == "StructuredColumns") {
    const auto fs =
        atlas::functionspace::StructuredColumns(xb[s].functionspace());
    if (fs.grid().name().compare(0, 1, "F") == 0 ||
        fs.grid().name().compare(0, 1, "O") == 0) {
      modelGrid = gauss;
    }
  } else if (xb[s].functionspace().type() == "NodeColumns") {
    const auto fs = atlas::functionspace::NodeColumns(xb[s].functionspace());
    if (fs.mesh().grid().name().compare(0, 2, "CS") == 0) {
      modelGrid = cubedsphere;
    }
  }
  if (modelGrid == unknown) {
    oops::Log::error() << "Only Gauss grid on StructuredColumns FunctionSpace"
                          " or cubed-sphere grid on NodeColumns FunctionSpace"
                          " are currently handled."
                       << std::endl;
    throw eckit::FunctionalityNotSupported(error_message, Here());
  }

  // Populate field "s"
  atlas::FieldSet modelFields;

  const bool lPopulateFields =
      (normfield(outerGeometryData.comm(), xb[s]) <= tolerance);

  if (lPopulateFields) {
    modelFields = populateFields(outerGeometryData.fieldSet(), xb);
  } else {
    modelFields.add(xb[s]);
  }
  modelFields[s].set_dirty();
  modelFields.haloExchange();

  // Interpolate to Gauss
  atlas::FieldSet gaussFields;
  if (modelGrid == cubedsphere) {
    interpolateCSToGauss(outerGeometryData, modelFields, gaussFields, normaliseInterpolation);
  } else if (modelGrid == gauss) {
    gaussFields = modelFields;
  }

  // Add Coriolis term
  gaussFields.add(createCoriolis(gaussFields[s]));

  return gaussFields;
}

// -----------------------------------------------------------------------------

atlas::FieldSet allocateSpectralVortDiv(
    const atlas::functionspace::Spectral& specfs, const std::size_t& levels) {
  atlas::FieldSet specfset;
  atlas::Field specvort = specfs.createField<double>(
      atlas::option::name("air_upward_absolute_vorticity") |
      atlas::option::levels(levels));
  atlas::Field specdiv = specfs.createField<double>(
      atlas::option::name("air_horizontal_divergence") |
      atlas::option::levels(levels));

  specfset.add(specvort);
  specfset.add(specdiv);

  return specfset;
}

// -----------------------------------------------------------------------------

oops::Variables removeOuterOnlyVar(const oops::Variables& vars) {
  oops::Variables innerVars(vars);
  innerVars -= innerVars["geostrophic_pressure_levels_minus_one"];
  return innerVars;
}

// -----------------------------------------------------------------------------

void applyRecipNtimesNplus1SpectralScaling(
    const oops::Variables& innerNames, const oops::Variables& outerNames,
    const atlas::functionspace::Spectral& specFS,
    const atlas::idx_t& totalWavenumber, atlas::FieldSet& fSet) {
  // copy fields that are not associated with innerNames
  atlas::FieldSet fsetTemp;
  for (atlas::Field& f : fSet) {
    if (!(innerNames.has(f.name()))) {
      fsetTemp.add(f);
    }
  }

  atlas::FieldSet fsetScaled;
  for (std::size_t var = 0; var < innerNames.size(); ++var) {
    atlas::Field scaledFld = fSet[innerNames[var].name()];
    const auto levels = scaledFld.shape(1);
    auto fldView = atlas::array::make_view<double, 2>(scaledFld);

    double earthRadius = atlas::util::Earth::radius();  // radius of earth;
    const double squaredEarthRadius = earthRadius * earthRadius;
    util::for_each_index(
        util::make_index_space_spectral_1d(scaledFld),
        [=](atlas::idx_t i, atlas::idx_t n1, atlas::idx_t m1) mutable {
          for (atlas::idx_t jl = 0; jl < levels; ++jl) {
            if (n1 != 0) {
              fldView(i, jl) *= squaredEarthRadius / (n1 * (n1 + 1));
            } else {
              fldView(i, jl) = 0.0;
            }
          }
        });

    scaledFld.rename(outerNames[var].name());
    fsetScaled.add(scaledFld);
  }

  fSet = fsetTemp;

  for (auto& f : fsetScaled) {
    fSet.add(f);
  }
}

}  //  namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussUVToGPWithSMV> makerGaussUVToGPWithSMV_(
    "gauss winds to geostrophic pressure with smv-interp");

GaussUVToGPWithSMV::GaussUVToGPWithSMV(
    const oops::GeometryData& outerGeometryData,
    const oops::Variables& outerVars, const eckit::Configuration& covarConf,
    const Parameters_& params, const oops::FieldSet3D& xb,
    const oops::FieldSet3D& fg)
    : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
      params_(params),
      innerVars_(removeOuterOnlyVar(
          getUnionOfInnerActiveAndOuterVars(params, outerVars))),
      activeOuterVars_(params.activeOuterVars(outerVars)),
      innerOnlyVars_(getInnerOnlyVars(params, outerVars)),
      gaussFunctionSpace_(outerGeometryData.functionSpace()),
      specFunctionSpace_(
          2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
      trans_(gaussFunctionSpace_, specFunctionSpace_),
      innerGeometryData_(atlas::FunctionSpace(gaussFunctionSpace_),
                         outerGeometryData.fieldSet(),
                         outerGeometryData.levelsAreTopDown(),
                         outerGeometryData.comm()),
      augmentedState_(createAugmentedState(outerGeometryData, xb.fieldSet(),
                               params.normalizeInverseInterpolation)) {
  oops::Log::trace() << classname() << "::GaussUVToGP starting" << std::endl;
  // read in "gaussian air density" state.
  oops::Log::trace() << classname() << "::GaussUVToGP done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGPWithSMV::multiply(oops::FieldSet3D& fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;
  // Allocate output fields if they are not already present, e.g when
  // randomizing.
  allocateMissingFields(fset, activeOuterVars_, activeOuterVars_,
                        innerGeometryData_.functionSpace());

  atlas::Field gp = fset["geostrophic_pressure_levels_minus_one"];

  atlas::Field rhsvec = allocateRHSVec(gaussFunctionSpace_, gp.shape(1));

  populateRHSVec(augmentedState_["dry_air_density_levels_minus_one"],
                 augmentedState_["coriolis"], fset.fieldSet(), rhsvec);

  atlas::FieldSet specfset =
      allocateSpectralVortDiv(specFunctionSpace_, rhsvec.shape(1));
  // calculate dir air_upward_absolute_vorticity and air_horizontal_divergence
  // spectrally
  trans_.dirtrans_wind2vordiv(rhsvec, specfset["air_upward_absolute_vorticity"],
                              specfset["air_horizontal_divergence"]);

  // apply inverse laplacian spectral scaling to spectral
  // air_horizontal_divergence
  const int N = specFunctionSpace_.truncation();
  applyRecipNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"air_horizontal_divergence"})),
      oops::Variables(std::vector<std::string>({"air_horizontal_divergence"})),
      specFunctionSpace_, N, specfset);

  // apply inverse spectral transform to
  trans_.invtrans(specfset["air_horizontal_divergence"], gp);

  gp.set_dirty();

  // Remove inner-only variables
  fset.removeFields(innerOnlyVars_);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGPWithSMV::multiplyAD(oops::FieldSet3D& fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  // Allocate inner-only variables
  checkFieldsAreNotAllocated(fset, innerOnlyVars_);
  allocateMissingFields(fset, innerOnlyVars_, innerOnlyVars_,
                        innerGeometryData_.functionSpace());

  atlas::FieldSet specfset = allocateSpectralVortDiv(
      specFunctionSpace_,
      fset["geostrophic_pressure_levels_minus_one"].shape(1));

  // apply inverse spectral transform to
  trans_.invtrans_adj(fset["geostrophic_pressure_levels_minus_one"],
                      specfset["air_horizontal_divergence"]);

  // apply inverse laplacian spectral scaling to spectral
  // air_horizontal_divergence
  const int N = specFunctionSpace_.truncation();
  applyRecipNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"air_horizontal_divergence"})),
      oops::Variables(std::vector<std::string>({"air_horizontal_divergence"})),
      specFunctionSpace_, N, specfset);

  auto vortView = atlas::array::make_view<double, 2>(
      specfset["air_upward_absolute_vorticity"]);
  vortView.assign(0.0);

  atlas::Field rhsvec =
      allocateRHSVec(gaussFunctionSpace_,
                     fset["geostrophic_pressure_levels_minus_one"].shape(1));

  trans_->dirtrans_wind2vordiv_adj(specfset["air_upward_absolute_vorticity"],
                                   specfset["air_horizontal_divergence"],
                                   rhsvec);

  populateRHSVecAdj(augmentedState_["dry_air_density_levels_minus_one"],
                    augmentedState_["coriolis"], fset.fieldSet(), rhsvec);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGPWithSMV::print(std::ostream& os) const { os << classname(); }

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
