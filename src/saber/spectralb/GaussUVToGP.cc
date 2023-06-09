/*
 * (C) Crown Copyright 2022-2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/GaussUVToGP.h"

#include <netcdf.h>

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

#include "mo/common_varchange.h"
#include "mo/eval_dry_air_density.h"
#include "mo/eval_virtual_potential_temperature.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/interpolation/AtlasInterpWrapper.h"

#include "saber/oops/SaberOuterBlockBase.h"

#define ERR(e) {ABORT(nc_strerror(e));}


namespace saber {
namespace spectralb {

namespace {

// -----------------------------------------------------------------------------

atlas::Field allocateRHSVec(const atlas::FunctionSpace & gaussFS,
                            const std::size_t & levels) {
  atlas::Field rhsvec = gaussFS.createField<double>(atlas::option::name("rhsvec") |
                                                    atlas::option::variables(2) |
                                                    atlas::option::levels(levels));
  auto rhsVecView = atlas::array::make_view<double, 3>(rhsvec);
  rhsVecView.assign(0.0);
  return rhsvec;
}

void populateRHSVec(const atlas::Field & rhoState,
                    const atlas::Field & coriolis,
                    const atlas::FieldSet & fset,
                    atlas::Field & rhsvec) {
  auto uView = atlas::array::make_view<const double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<const double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto coriolisView = atlas::array::make_view<const double, 2>(coriolis);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  for (atlas::idx_t jn = 0; jn < rhsvecView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < rhsvecView.shape()[1]; ++jl) {
      rhsvecView(jn, jl, atlas::LON) =
        - coriolisView(jn, 0) * rhoStateView(jn, jl) * vView(jn, jl);
      rhsvecView(jn, jl, atlas::LAT) =
        coriolisView(jn, 0) * rhoStateView(jn, jl) * uView(jn, jl);
    }
  }
}

void populateRHSVecAdj(const atlas::Field & rhoState,
                       const atlas::Field & coriolis,
                       atlas::FieldSet & fset,
                       atlas::Field & rhsvec) {
  auto uView = atlas::array::make_view<double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto coriolisView = atlas::array::make_view<const double, 2>(coriolis);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  for (atlas::idx_t jn = 0; jn < rhsvecView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < rhsvecView.shape()[1]; ++jl) {
      vView(jn, jl) += - coriolisView(jn, 0) * rhoStateView(jn, jl)
                       * rhsvecView(jn, jl, atlas::LON);
      uView(jn, jl) += coriolisView(jn, 0) * rhoStateView(jn, jl)
                       * rhsvecView(jn, jl, atlas::LAT);
    }
  }
  rhsvecView.assign(0.0);
}

double normfield(const eckit::mpi::Comm & comm, const atlas::Field& fld) {
  // to do - maybe not include halos?
  auto view = atlas::array::make_view<double, 2>(fld);

  std::size_t i(0);
  double zz(0.0);
  for (atlas::idx_t jn = 0; jn < fld.shape(0); ++jn) {
    for (atlas::idx_t jl = 0; jl < fld.shape(1); ++jl, ++i) {
      zz += view(jn, jl) * view(jn, jl);
    }
  }
  comm.allReduceInPlace(i, eckit::mpi::sum());
  comm.allReduceInPlace(zz, eckit::mpi::sum());

  return std::sqrt(zz)/static_cast<double>(i);
}


void populateFields(const atlas::FieldSet & geomfields,
                    const atlas::FieldSet & statefields,
                    atlas::FieldSet & outputfields) {
  atlas::FieldSet tempfields;

  // Need to setup derived state fields that we need (done on model grid).
  std::vector<std::string> requiredStateVariables{
    "air_pressure_levels_minus_one",
    "dry_air_density_levels_minus_one",
    "height",
    "height_levels",
    "m_ci",
    "m_cl",
    "m_r",
    "m_v",
    "m_t",
    "potential_temperature",
    "specific_humidity",
    "virtual_potential_temperature"};

  std::vector<std::string> requiredGeometryVariables{"height_levels",
                                                     "height"};

  std::vector<std::string> outputVariables{"dry_air_density_levels_minus_one"};

  // Check that they are allocated (i.e. exist in the state fieldset)
  for (auto & s : requiredStateVariables) {
    if (!statefields.has(s)) {
      oops::Log::error() << "::gaussuvtogp::populateFields variable " << s <<
                            "is not part of state object." << std::endl;
    }
  }

  tempfields.clear();
  outputfields.clear();
  for (const auto & s : requiredStateVariables) {
    tempfields.add(statefields[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    tempfields.add(geomfields[s]);
  }

  mo::evalTotalMassMoistAir(tempfields);
  mo::evalSpecificHumidity(tempfields);
  mo::eval_virtual_potential_temperature_nl(tempfields);
  mo::eval_dry_air_density_nl(tempfields);

  std::string s("dry_air_density_levels_minus_one");
  outputfields.add(tempfields[s]);
}

void interpolateCSToGauss(const std::string & modelGridName,
                          const oops::GeometryData & outerGeometryData,
                          const atlas::FieldSet &  csfields,
                          atlas::FieldSet & gfields) {
  std::string s("dry_air_density_levels_minus_one");

  const auto srcFunctionspace =
      atlas::functionspace::NodeColumns(csfields[s].functionspace());

  atlas::CubedSphereGrid grid(modelGridName);
  const auto meshConfig = atlas::util::Config("partitioner", "cubedsphere")
    | atlas::util::Config("halo", 0);
  const auto meshGen = atlas::MeshGenerator("cubedsphere_dual", meshConfig);
  auto csmesh = atlas::Mesh(meshGen.generate(grid));

  const auto partitioner = atlas::grid::MatchingPartitioner(srcFunctionspace.mesh(),
                           atlas::util::Config("type", "cubedsphere"));

  const auto targetGrid =
    atlas::functionspace::StructuredColumns(outerGeometryData.functionSpace()).grid();
  const auto targetMesh = atlas::MeshGenerator("structured").generate(targetGrid,
                                                                      partitioner);
  const auto targetFunctionspace = atlas::functionspace::NodeColumns(targetMesh);

  // Set up interpolation object.
  const auto scheme = atlas::util::Config("type", "cubedsphere-bilinear") |
                      atlas::util::Config("adjoint", false);
  const auto interp = atlas::Interpolation(scheme, csfields[0].functionspace(),
                                           targetFunctionspace);

  auto targetField = targetFunctionspace.createField<double>
    (atlas::option::name(s) |
     atlas::option::levels(csfields[s].levels()) |
                           atlas::option::halo(1));

  interp.execute(csfields[s], targetField);

  auto atlaswrap = ::saber::interpolation::AtlasInterpWrapper(
    partitioner, targetFunctionspace,
    targetGrid, outerGeometryData.functionSpace());

  atlas::Field gaussField =
    atlas::functionspace::StructuredColumns(outerGeometryData.functionSpace()).createField<double>(
    (atlas::option::name(s) | atlas::option::levels(csfields[s].levels()) |
     atlas::option::halo(1)));

  atlaswrap.execute(targetField, gaussField);

  gaussField.haloExchange();

  gfields.add(gaussField);
}

atlas::Field createCoriolis(const atlas::Field & scstate) {
  // create coriolis parameter.
  // need to get the latitude at each grid point from the mesh.
  static constexpr double twoOmega = 2.0 * 7.292116E-5;
  atlas::Field coriolis = scstate.functionspace().createField<double>(
    atlas::option::name("coriolis") | atlas::option::levels(1));

  auto coriolisView = atlas::array::make_view<double, 2>(coriolis);
  auto sc = atlas::functionspace::StructuredColumns(scstate.functionspace());
  atlas::idx_t jn(0);
  for (atlas::idx_t j = sc.j_begin(); j < sc.j_end(); ++j) {
    for (atlas::idx_t i = sc.i_begin(j); i < sc.i_end(j); ++i) {
      jn = sc.index(i, j);
      coriolisView(jn, 0) = twoOmega * sin(atlas::util::Constants::degreesToRadians() *
                                           sc.grid().lonlat(i, j).lat());
    }
  }
  coriolis.haloExchange();

  return coriolis;
}

atlas::FieldSet createAugmentedState(const std::string & modelGridName,
                                     const std::string & gaussStateName,
                                     const oops::GeometryData & outerGeometryData,
                                     const atlas::FieldSet & xb) {
  atlas::FieldSet tempfields;
  atlas::FieldSet gfields;

  std::string s("dry_air_density_levels_minus_one");
  double tolerance(1.0e-6);

  // error trap that xb has field "s" even if it is not populated.
  if (!xb.has(s)) {
    oops::Log::error() <<  "expect " << s
                       << "to be allocated in xb fieldset"
                       << std::endl;
    throw std::runtime_error("ERROR - gaussuvtogp saber block: failed");
  }

  if ( (xb[s].functionspace().type() == "StructuredColumns") &&
       (atlas::functionspace::StructuredColumns(
          xb[s].functionspace()).grid().name().compare(0, 1, std::string("F")) == 0) ) {
    if (normfield(outerGeometryData.comm(), xb[s]) <= tolerance) {
      populateFields(outerGeometryData.fieldSet(), xb, gfields);
    } else {
      gfields.add(xb[s]);
    }
    gfields.add(createCoriolis(xb[s]));
    return gfields;
  }

  if (gaussStateName.size() > 0)  {
    atlas::FieldSet globalData;

    atlas::Field field = outerGeometryData.functionSpace().createField<double>(
      atlas::option::name(s) | atlas::option::levels(xb[s].levels()) |
      atlas::option::global());
    globalData.add(field);

    // read in gauss state on a single PE.
    if (outerGeometryData.comm().rank() == 0) {
      atlas::StructuredGrid grid =
        atlas::functionspace::StructuredColumns(outerGeometryData.functionSpace()).grid();

      // Get sizes
      atlas::idx_t nx = grid.nxmax();
      atlas::idx_t ny = grid.ny();

      // NetCDF IDs
      int ncid, retval, var_id;

      // NetCDF file path
      std::string ncfilepath = gaussStateName;
      oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

      // Open NetCDF file
      if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

      oops::Log::info() << "Info     : Reading file done: " << retval << std::endl;

      // Get variable
      if ((retval = nc_inq_varid(ncid, s.c_str(), &var_id))) ERR(retval);

      // Read data
      double zvar[xb[s].levels()][ny][nx];

      if ((retval = nc_get_var_double(ncid, var_id, &zvar[0][0][0]))) ERR(retval);
      // Copy data

      auto varView = atlas::array::make_view<double, 2>(globalData[s]);
      for (atlas::idx_t k = 0; k < xb[s].levels(); ++k) {
        for (atlas::idx_t j = 0; j < ny; ++j) {
          for (atlas::idx_t i = 0; i < grid.nx(ny-1-j); ++i) {
            atlas::gidx_t gidx = grid.index(i, ny-1-j);
            varView(gidx, k) = zvar[k][j][i];
          }
        }
      }
      // Close file
      if ((retval = nc_close(ncid))) ERR(retval);
    }

    // redistribute field across PEs.
    atlas::functionspace::StructuredColumns fs(outerGeometryData.functionSpace());

    atlas::Field fieldLoc =
      outerGeometryData.functionSpace().createField<double>(atlas::option::name(s)
      | atlas::option::levels(xb[s].levels()) | atlas::option::halo(1));
    gfields.add(fieldLoc);
    fs.scatter(globalData, gfields);

    // create Coriolis
    gfields.add(createCoriolis(gfields[s]));

    // halo exchange
    gfields.haloExchange();

    return gfields;
  }

  // We don't have a gauss state to read in, so we interpolate
  // Unfortunately does not work all the way to StructureColumns
  // ie only valid for a single2 1PE (for now)

  if (normfield(outerGeometryData.comm(), xb[s]) > tolerance) {
    // just link to existing data.
    tempfields.add(xb[s]);
  } else {
    // gauss field exists but is not populated.
    populateFields(outerGeometryData.fieldSet(), xb, tempfields);
  }

  if ((xb[s].functionspace().type() == "NodeColumns") &&
      (modelGridName.compare(0, 2, std::string{"CS"}) == 0)) {
    interpolateCSToGauss(modelGridName, outerGeometryData, tempfields, gfields);
  } else {
    gfields.add(tempfields[s]);
  }
  // create Coriolis
  gfields.add(createCoriolis(gfields[s]));

  return gfields;
}


// -----------------------------------------------------------------------------

atlas::FieldSet allocateSpectralVortDiv(
    const atlas::functionspace::Spectral & specfs,
    const std::size_t & levels) {
  atlas::FieldSet specfset;
  atlas::Field specvort = specfs.createField<double>(
    atlas::option::name("vorticity") | atlas::option::levels(levels));
  atlas::Field specdiv = specfs.createField<double>(
    atlas::option::name("divergence") | atlas::option::levels(levels));

  specfset.add(specvort);
  specfset.add(specdiv);

  return specfset;
}


// -----------------------------------------------------------------------------

oops::Variables createInnerVars(const oops::Variables & outerVars) {
  oops::Variables innerVars;

  for (auto & var : outerVars.variables()) {
    if (var.compare("geostrophic_pressure_levels_minus_one") == 0) {
    } else {
      innerVars.push_back(var);
    }
  }
  return innerVars;
}

// -----------------------------------------------------------------------------

void applyRecipNtimesNplus1SpectralScaling(const oops::Variables & innerNames,
                                           const oops::Variables & outerNames,
                                           const atlas::functionspace::Spectral & specFS,
                                           const atlas::idx_t & totalWavenumber,
                                           atlas::FieldSet & fSet) {
  const auto zonal_wavenumbers = specFS.zonal_wavenumbers();
  const int nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // copy fields that are not associated with innerNames
  atlas::FieldSet fsetTemp;
  for (atlas::Field & f : fSet) {
    if (!(innerNames.has(f.name()))) {
      fsetTemp.add(f);
    }
  }

  atlas::FieldSet fsetScaled;
  for (std::size_t var = 0; var < innerNames.variables().size(); ++var) {
    atlas::Field scaledFld = fSet[innerNames[var]];
    auto fldView = atlas::array::make_view<double, 2>(scaledFld);

    double earthRadius = atlas::util::Earth::radius();  // radius of earth;
    const double squaredEarthRadius = earthRadius * earthRadius;
    int i(0);
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(totalWavenumber); ++n1) {
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (atlas::idx_t jl = 0; jl < fSet[innerNames[var]].levels(); ++jl) {
            if (n1 != 0) {
              fldView(i, jl) *= squaredEarthRadius / (n1 * (n1 + 1));
            } else {
              fldView(i, jl) = 0.0;
            }
          }
        }
      }
    }
    scaledFld.rename(outerNames[var]);
    fsetScaled.add(scaledFld);
  }

  fSet = fsetTemp;

  for (auto & f : fsetScaled) {
    fSet.add(f);
  }
}

}  //  namespace

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<GaussUVToGP>
   makerGaussUVToGP_("gauss winds to geostrophic pressure");


GaussUVToGP::GaussUVToGP(const oops::GeometryData & outerGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const oops::Variables & outerVars,
               const eckit::Configuration & covarConf,
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    params_(params),
    innerVars_(createInnerVars(outerVars)),
    outerVars_(outerVars),
    activeVariableSizes_(activeVariableSizes),
    modelGridName_(params_.modelGridName.value().get_value_or("")),
    gaussStateName_(params_.gaussState.value().get_value_or("")),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    innerGeometryData_(atlas::FunctionSpace(gaussFunctionSpace_),
                       outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(), outerGeometryData.comm()),
    augmentedState_(createAugmentedState(modelGridName_, gaussStateName_,
                                         outerGeometryData, xb))
{
  oops::Log::trace() << classname() << "::GaussUVToGP starting" << std::endl;
  // read in "gaussian air density" state.
  oops::Log::trace() << classname() << "::GaussUVToGP done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << std::endl;

  atlas::Field gp = gaussFunctionSpace_.createField<double>(
    atlas::option::name("geostrophic_pressure_levels_minus_one") |
    atlas::option::levels(fset["eastward_wind"].levels()));

  atlas::Field rhsvec = allocateRHSVec(gaussFunctionSpace_, gp.levels());

  populateRHSVec(augmentedState_["dry_air_density_levels_minus_one"],
                 augmentedState_["coriolis"], fset, rhsvec);

  atlas::FieldSet specfset = allocateSpectralVortDiv(specFunctionSpace_, rhsvec.levels());
  // calculate dir vorticity and divergence spectrally
  trans_.dirtrans_wind2vordiv(rhsvec, specfset["vorticity"], specfset["divergence"]);

  // apply inverse laplacian spectral scaling to spectral divergence
  const int N = specFunctionSpace_.truncation();
  applyRecipNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"divergence"})),
      oops::Variables(std::vector<std::string>({"divergence"})),
      specFunctionSpace_, N, specfset);

  // apply inverse spectral transform to
  trans_.invtrans(specfset["divergence"], gp);

  gp.haloExchange();

  if (fset.has("geostrophic_pressure_levels_minus_one")) {
    auto gpView =
      atlas::array::make_view<double, 2>(
        fset["geostrophic_pressure_levels_minus_one"]);
    auto gpViewKeep =
      atlas::array::make_view<const double, 2>(gp);
    gpView.assign(gpViewKeep);
  } else {
    fset.add(gp);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  atlas::FieldSet specfset =
      allocateSpectralVortDiv(specFunctionSpace_,
                              fset["geostrophic_pressure_levels_minus_one"].levels());

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  oops::Variables
    activeVars(std::vector<std::string>({"eastward_wind", "northward_wind",
                                         "geostrophic_pressure_levels_minus_one"}));

  for (auto & s : fset.field_names()) {
    if (!activeVars.has(s)) {
      newFields.add(fset[s]);
    }
  }

  fset["geostrophic_pressure_levels_minus_one"].adjointHaloExchange();

  // apply inverse spectral transform to
  trans_.invtrans_adj(fset["geostrophic_pressure_levels_minus_one"], specfset["divergence"]);

  // apply inverse laplacian spectral scaling to spectral divergence
  const int N = specFunctionSpace_.truncation();
  applyRecipNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"divergence"})),
      oops::Variables(std::vector<std::string>({"divergence"})),
      specFunctionSpace_, N, specfset);

  auto vortView = atlas::array::make_view<double, 2>(specfset["vorticity"]);
  vortView.assign(0.0);

  atlas::Field rhsvec = allocateRHSVec(gaussFunctionSpace_,
                                       fset["geostrophic_pressure_levels_minus_one"].levels());

  trans_->dirtrans_wind2vordiv_adj(specfset["vorticity"], specfset["divergence"], rhsvec);

  populateRHSVecAdj(augmentedState_["dry_air_density_levels_minus_one"],
                    augmentedState_["coriolis"],
                    fset, rhsvec);

  newFields.add(fset["eastward_wind"]);
  newFields.add(fset["northward_wind"]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
