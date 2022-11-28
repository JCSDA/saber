/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralToGaussUV.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Earth.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace spectralb {

namespace {

// -----------------------------------------------------------------------------

atlas::Field allocateGaussUVField(const atlas::FunctionSpace & gaussFS,
                                  const oops::Variables & innerVariables,
                                  const std::vector<std::size_t> & activeVariableSizes) {
  std::array<size_t, 2> indx{{0, 0}};
  std::array<size_t, 2> levels{{0, 0}};
  if (innerVariables.has("vorticity") && innerVariables.has("divergence")) {
    indx[0] = innerVariables.find("vorticity");
    indx[1] = innerVariables.find("divergence");
  } else if (innerVariables.has("streamfunction") && innerVariables.has("velocity_potential")) {
    indx[0] = innerVariables.find("streamfunction");
    indx[1] = innerVariables.find("velocity_potential");
  } else {
    // error trap
    oops::Log::error() << "ERROR - either vorticity and divergence "
                       << "or streamfunction and velocity_potential "
                       << "not present " << std::endl;
    throw std::runtime_error("inner fields mis-specified");
  }
  // check that levels in indx[0] and indx[1] the same.

  levels[0] = activeVariableSizes[indx[0]];
  levels[1] = activeVariableSizes[indx[1]];

  if (levels[0] != levels[1]) {
    oops::Log::error() << "ERROR - the number of model levels in "
                       << "vorticity and divergence or "
                       << "streamfunction and velocity potential "
                       << indx[0] << " "
                       << indx[1] << " "
                       << levels[0] << " "
                       << levels[1]
                       << std::endl;
    throw std::runtime_error("vertical levels are inconsistent");
  }

  atlas::idx_t modellevels = static_cast<atlas::idx_t>(levels[0]);
  auto sc = atlas::functionspace::StructuredColumns(gaussFS);
  atlas::Field uvgp = sc.createField<double>(atlas::option::name("uv_gp") |
                                             atlas::option::variables(2) |
                                             atlas::option::levels(modellevels));
  return uvgp;
}

// -----------------------------------------------------------------------------

atlas::FieldSet allocateSpectralVortDiv(
    const atlas::functionspace::Spectral & specfs,
    const oops::Variables & activeVariables,
    const std::vector<std::size_t> & activeVariableSizes) {
  std::array<size_t, 2> indx;
  std::array<size_t, 2> levels;


  if (activeVariables.has("vorticity") && activeVariables.has("divergence")) {
    indx[0] = activeVariables.find("vorticity");
    indx[1] = activeVariables.find("divergence");
  } else if (activeVariables.has("streamfunction") && activeVariables.has("velocity_potential")) {
    indx[0] = activeVariables.find("streamfunction");
    indx[1] = activeVariables.find("velocity_potential");
  } else {
    // error trap
    oops::Log::error() << "ERROR - either vorticity and divergence "
                       << "or streamfunction and velocity potential "
                       << "not present " << std::endl;
    throw std::runtime_error("inner fields mis-specified");
  }
  // check that levels in indx[0] and indx[1] the same.
  levels[0] = activeVariableSizes[indx[0]];
  levels[1] = activeVariableSizes[indx[1]];


  atlas::FieldSet specfset;
  atlas::Field specvort = specfs.createField<double>(
    atlas::option::name("vorticity") | atlas::option::levels(levels[0]));
  atlas::Field specdiv = specfs.createField<double>(
    atlas::option::name("divergence") | atlas::option::levels(levels[1]));

  specfset.add(specvort);
  specfset.add(specdiv);

  return specfset;
}

// -----------------------------------------------------------------------------

atlas::FieldSet convertUVToFieldSet(const atlas::Field & uvField) {
  atlas::FieldSet uvfset;
  // It might make sense to either have
  //     both uvField and uvFieldSet having a halo.
  //     have uvField without halo and uvFieldSet with
  // I am going with the latter.

  atlas::Field u = uvField.functionspace().createField<double>
      (atlas::option::name("eastward_wind") |
       atlas::option::levels(uvField.levels()));

  atlas::Field v = uvField.functionspace().createField<double>
      (atlas::option::name("northward_wind") |
       atlas::option::levels(uvField.levels()));

  auto uView = atlas::array::make_view<double, 2>(u);
  auto vView = atlas::array::make_view<double, 2>(v);
  auto uvView = atlas::array::make_view<double, 3>(uvField);

  for (atlas::idx_t jn = 0; jn < uvView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < uvView.shape()[1]; ++jl) {
      uView(jn, jl) = uvView(jn, jl, 0);
      vView(jn, jl) = uvView(jn, jl, 1);
    }
  }
  u.haloExchange();
  v.haloExchange();

  uvfset.add(u);
  uvfset.add(v);

  return uvfset;
}

// -----------------------------------------------------------------------------

atlas::Field convertUVToFieldSetAD(const atlas::FieldSet & fset) {
  atlas::idx_t modellevels = static_cast<atlas::idx_t>(fset["eastward_wind"].levels());
  auto sc = atlas::functionspace::StructuredColumns(fset["eastward_wind"].functionspace());

  atlas::Field uvgp = sc.createField<double>(atlas::option::name("uv_gp") |
                                             atlas::option::variables(2) |
                                             atlas::option::levels(modellevels));


  fset["eastward_wind"].adjointHaloExchange();
  fset["northward_wind"].adjointHaloExchange();

  auto uView = atlas::array::make_view<double, 2>(fset["eastward_wind"]);
  auto vView = atlas::array::make_view<double, 2>(fset["northward_wind"]);
  auto uvView = atlas::array::make_view<double, 3>(uvgp);

  for (atlas::idx_t jn = 0; jn < uvView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < uvView.shape()[1]; ++jl) {
      uvView(jn, jl, 0) = uView(jn, jl);
      uvView(jn, jl, 1) = vView(jn, jl);
    }
  }

  return uvgp;
}

// -----------------------------------------------------------------------------

// Create inner variables from outer variables, excluding winds, and adding either
// (vorticity and divergence) or (stream function and velocity potential), depending on
// the value of params.useStreamFunctionVelocityPotential.
oops::Variables createInnerVars(const SpectralToGaussUVParameters & params,
                                const oops::Variables & outerVars) {
  oops::Variables innerVars;

  for (auto & var : outerVars.variables()) {
    if (var.compare("eastward_wind") == 0 || var.compare("northward_wind") == 0) {
    } else {
      innerVars.push_back(var);
    }
  }

  if (params.useStreamFunctionVelocityPotential.value()) {
    innerVars.push_back("streamfunction");
    innerVars.push_back("velocity_potential");
  } else {
    innerVars.push_back("vorticity");
    innerVars.push_back("divergence");
  }
  return innerVars;
}

// -----------------------------------------------------------------------------

void applyNtimesNplus1SpectralScaling(const oops::Variables & innerNames,
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
            fldView(i, jl) =  n1 * (n1 + 1) *  fldView(i, jl) / (a * a);
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

static SaberOuterBlockMaker<SpectralToGaussUV> makerSpectralToGaussUV_("spectral to gauss winds");

// Build inner functionspace from outer functionspace
// It is the outer functionspace that is in the argument.
// So we need to create spectral functionspace !!!

// We need to write separate functions that allow

// Rule:
// Expect each saber block in its multiply method as a final step to apply a
// halo exchange at the end of the method if it needs it.
// and at the beginning of the adjoint of the saber block to apply the adjoint
// halo exchange.

// "active variables" - now required in yaml
// "inner variables" - optional in yaml - inner for the multiply
//                   - sum active and passive variables
SpectralToGaussUV::SpectralToGaussUV(const oops::GeometryData & outerGeometryData,
               const std::vector<size_t> & activeVariableSizes,
               const oops::Variables & outerVars,
               const Parameters_ & params,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : params_(params),
    innerVars_(createInnerVars(params_, outerVars)),
    outerVars_(outerVars),
    activeVariableSizes_(activeVariableSizes),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    innerGeometryData_((params.useInnerGaussianFunctionSpace.value() ?
                        atlas::FunctionSpace(gaussFunctionSpace_) :
                        atlas::FunctionSpace(specFunctionSpace_)),
                        outerGeometryData.fieldSet(),
                        outerGeometryData.levelsAreTopDown(), outerGeometryData.comm())
{
  oops::Log::trace() << classname() << "::SpectralToGaussUV starting" << std::endl;
  oops::Log::trace() << classname() << "::SpectralToGaussUV done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGaussUV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting " << fset.field_names() << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  for (auto & s : fset.field_names()) {
     if (outerVars_.has(s)) {
       newFields.add(fset[s]);
     }
  }

  // check variables and if streamfunction, velocity potential scale
  // by n(n+1) / a   and rename fields to vorticity_spectral_2D
  // divergence_spectral_2D
  const int N = specFunctionSpace_.truncation();

  if (innerVars_.has("streamfunction") && innerVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"streamfunction", "velocity_potential"})),
      oops::Variables(std::vector<std::string>({"vorticity", "divergence"})),
      specFunctionSpace_, N, fset);
  }

  atlas::Field uvgp = allocateGaussUVField(gaussFunctionSpace_,
                                           innerVars_, activeVariableSizes_);

  // transform to gaussian grid
  trans_.invtrans_vordiv2wind(fset["vorticity"], fset["divergence"], uvgp);

  atlas::FieldSet uvfset = convertUVToFieldSet(uvgp);

  newFields.add(uvfset["eastward_wind"]);
  newFields.add(uvfset["northward_wind"]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGaussUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  for (auto & s : fset.field_names()) {
     if (innerVars_.has(s)) {
       newFields.add(fset[s]);
     }
  }

  atlas::Field uvgp = convertUVToFieldSetAD(fset);
  atlas::FieldSet specfset = allocateSpectralVortDiv(specFunctionSpace_,
                                                     innerVars_, activeVariableSizes_);

  trans_.invtrans_vordiv2wind_adj(uvgp, specfset["vorticity"], specfset["divergence"]);

  const int N = specFunctionSpace_.truncation();

  if (innerVars_.has("streamfunction") && innerVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"vorticity", "divergence"})),
      oops::Variables(std::vector<std::string>({"streamfunction", "velocity_potential"})),
      specFunctionSpace_, N, specfset);
  }

  newFields.add(specfset[0]);
  newFields.add(specfset[1]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGaussUV::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGaussUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
