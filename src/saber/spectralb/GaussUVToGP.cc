/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/GaussUVToGP.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Earth.h"

#include "mo/control2analysis_varchange.h"
#include "mo/model2geovals_varchange.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

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
                    const atlas::FieldSet & fset,
                    atlas::Field & rhsvec) {
  auto uView = atlas::array::make_view<const double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<const double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  for (atlas::idx_t jn = 0; jn < rhsvecView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < rhsvecView.shape()[1]; ++jl) {
      rhsvecView(jn, jl, atlas::LON) = - rhoStateView(jn, jl) * vView(jn, jl);
      rhsvecView(jn, jl, atlas::LAT) = rhoStateView(jn, jl) * uView(jn, jl);
    }
  }
}

void populateRHSVecAdj(const atlas::Field & rhoState,
                       atlas::FieldSet & fset,
                       atlas::Field & rhsvec) {
  auto uView = atlas::array::make_view<double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<double, 2>(fset["northward_wind"]);
  auto rhoStateView = atlas::array::make_view<const double, 2>(rhoState);
  auto rhsvecView = atlas::array::make_view<double, 3>(rhsvec);

  for (atlas::idx_t jn = 0; jn < rhsvecView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < rhsvecView.shape()[1]; ++jl) {
      vView(jn, jl) += - rhoStateView(jn, jl) * rhsvecView(jn, jl, atlas::LON);
      uView(jn, jl) += rhoStateView(jn, jl) * rhsvecView(jn, jl, atlas::LAT);
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
  std::vector<std::string> requiredStateVariables{ "exner_levels_minus_one",
                                                   "potential_temperature",
                                                   "exner",
                                                   "air_pressure_levels_minus_one",
                                                   "air_temperature",
                                                   "dry_air_density_levels_minus_one"};

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
  for (const auto & s : requiredStateVariables) {
    tempfields.add(statefields[s]);
  }

  for (const auto & s : requiredGeometryVariables) {
    tempfields.add(geomfields[s]);
  }

  mo::evalAirTemperature(tempfields);
  mo::evalDryAirDensity(tempfields);

  std::string s("dry_air_density_levels_minus_one");
  outputfields.add(tempfields[s]);
}

void interpolateCSToGauss(const atlas::Grid & modelGrid,
                          const oops::GeometryData & outerGeometryData,
                          const atlas::FieldSet &  csfields,
                          atlas::FieldSet & gfields) {
  const auto srcFunctionspace = atlas::functionspace::NodeColumns(csfields[0].functionspace());

  const auto partitioner = atlas::grid::MatchingPartitioner(srcFunctionspace,
                           atlas::util::Config("type", "cubedsphere"));

  const auto targetGrid =
      atlas::functionspace::StructuredColumns(outerGeometryData.functionSpace()).grid();
  const auto targetMesh = atlas::MeshGenerator("structured").generate(targetGrid, partitioner);
  const auto targetFunctionspace = atlas::functionspace::NodeColumns(targetMesh);

  // Set up interpolation object.
  const auto scheme = atlas::util::Config("type", "cubedsphere-bilinear") |
                      atlas::util::Config("adjoint", true);
  const auto interp = atlas::Interpolation(scheme, csfields[0].functionspace(),
                                           targetFunctionspace);
  std::string s("dry_air_density_levels_minus_one");
  auto targetField = targetFunctionspace.createField<double>
      (atlas::option::name(s) |
       atlas::option::levels(csfields[s].levels()) |
                             atlas::option::halo(1));
  gfields.add(targetField);
}


atlas::FieldSet createAugmentedState(const oops::GeometryData & outerGeometryData,
                                     const atlas::FieldSet & xb) {
  // below will work if xb is on the gauss mesh.
  // maybe error trap to check functionspace?
  // maybe alternative is to read in GaussField from file.

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

  if (xb[0].functionspace().type() == "StructuredColumns") {
    const auto srcFunctionspace =
      atlas::functionspace::StructuredColumns(xb[0].functionspace());
    if (srcFunctionspace.grid().name().compare(0, 1, std::string{"F"}) == 0) {
      if (normfield(outerGeometryData.comm(), xb[s]) > tolerance) {
        // just link to existing data.
        gfields.add(xb[s]);
      } else {
        // gauss field exists but is not populated.
        populateFields(outerGeometryData.fieldSet(), xb, gfields);
      }
    }
  } else if (xb[0].functionspace().type() == "NodeColumns") {
    const auto srcFunctionspace = atlas::functionspace::NodeColumns(xb[0].functionspace());
    if (srcFunctionspace.mesh().grid().name().compare(0, 2, std::string{"CS"}) == 0) {
      if (normfield(outerGeometryData.comm(), xb[s]) <= tolerance) {
          populateFields(outerGeometryData.fieldSet(), xb, tempfields);
      }
      interpolateCSToGauss(outerGeometryData, tempfields, gfields);
    } else {
      oops::Log::error() <<  "non cubed-sphere NodeColumns not supported "
                         << std::endl;
      throw std::runtime_error("ERROR - gaussuvtogp saber block: failed");
    }
  } else {
    oops::Log::error() <<  "functionspace type "
                       << xb[0].functionspace().type()
                       << "not currently supported"
                       << std::endl;
    throw std::runtime_error("ERROR - gaussuvtogp saber block: failed");
  }

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
    if (var.compare("geostrophic_pressure") == 0) {
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

    double a = atlas::util::Earth::radius();  // radius of earth;
    int i(0);
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(totalWavenumber); ++n1) {
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (atlas::idx_t jl = 0; jl < fSet[innerNames[var]].levels(); ++jl) {
            if (n1 != 0) {
              fldView(i, jl) = a * a *  fldView(i, jl) / (n1 * (n1 + 1));
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
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : params_(params),
    innerVars_(createInnerVars(outerVars)),
    outerVars_(outerVars),
    activeVariableSizes_(activeVariableSizes),
    modelgrid_(params_.modelGridName),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    innerGeometryData_(atlas::FunctionSpace(gaussFunctionSpace_),
                       outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(), outerGeometryData.comm()),
    augmentedState_(createAugmentedState(modelgrid_, outerGeometryData, xb))
{
  oops::Log::trace() << classname() << "::GaussUVToGP starting" << std::endl;
  // read in "gaussian air density" state.
  oops::Log::trace() << classname() << "::GaussUVToGP done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << fset.field_names() << std::endl;

  atlas::Field gp = gaussFunctionSpace_.createField<double>(
    atlas::option::name("geostrophic_pressure") |
    atlas::option::levels(fset["eastward_wind"].levels()) |
    atlas::option::halo(1));

  // calculate 2D vector vec
  atlas::Field rhsvec = allocateRHSVec(gaussFunctionSpace_, gp.levels());

  populateRHSVec(augmentedState_["dry_air_density_levels_minus_one"], fset, rhsvec);

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

  fset.add(gp);

  oops::Log::trace() << classname() << "::multiply done"
                     << fset.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << fset.field_names() <<std::endl;

  atlas::FieldSet specfset = allocateSpectralVortDiv(specFunctionSpace_,
                                                     fset["geostrophic_pressure"].levels());

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  for (auto & s : fset.field_names()) {
     if (!innerVars_.has(s) && (s.compare("geostrophic_pressure") != 0)) {
       newFields.add(fset[s]);
     }
  }
  fset["geostrophic_pressure"].adjointHaloExchange();

  // apply inverse spectral transform to
  trans_.invtrans_adj(fset["geostrophic_pressure"], specfset["divergence"]);

  // apply inverse laplacian spectral scaling to spectral divergence
  const int N = specFunctionSpace_.truncation();
  applyRecipNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"divergence"})),
      oops::Variables(std::vector<std::string>({"divergence"})),
      specFunctionSpace_, N, specfset);

  auto vortView = atlas::array::make_view<double, 2>(specfset["vorticity"]);
  vortView.assign(0.0);

  atlas::Field rhsvec = allocateRHSVec(gaussFunctionSpace_,
                                       fset["geostrophic_pressure"].levels());

  populateRHSVec(augmentedState_["dry_air_density_levels_minus_one"], fset, rhsvec);

  trans_->dirtrans_wind2vordiv_adj(specfset["vorticity"], specfset["divergence"], rhsvec);

  populateRHSVecAdj(augmentedState_["dry_air_density_levels_minus_one"], fset, rhsvec);

  newFields.add(fset["eastward_wind"]);
  newFields.add(fset["northward_wind"]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << fset.field_names() << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void GaussUVToGP::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
