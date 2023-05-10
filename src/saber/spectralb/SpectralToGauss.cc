/*
 * (C) Copyright 2022- UCAR
 * (C) Crown Copyright 2022- Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SpectralToGauss.h"

#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

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

  const atlas::idx_t modellevels = static_cast<atlas::idx_t>(levels[0]);
  const auto sc = atlas::functionspace::StructuredColumns(gaussFS);
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
  } else if (activeVariables.has("streamfunction") &&
             activeVariables.has("velocity_potential")) {
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
  const auto uvView = atlas::array::make_view<double, 3>(uvField);

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
  const atlas::Field & uField = fset["eastward_wind"];
  const atlas::Field & vField = fset["northward_wind"];

  atlas::idx_t modellevels = static_cast<atlas::idx_t>(uField.levels());
  const auto sc = atlas::functionspace::StructuredColumns(uField.functionspace());

  atlas::Field uvgp = sc.createField<double>(atlas::option::name("uv_gp") |
                                             atlas::option::variables(2) |
                                             atlas::option::levels(modellevels));

  uField.adjointHaloExchange();
  vField.adjointHaloExchange();

  const auto uView = atlas::array::make_view<double, 2>(uField);
  const auto vView = atlas::array::make_view<double, 2>(vField);
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

atlas::Field convertFieldSetToUV(const atlas::FieldSet & fset) {
  const atlas::Field & uField = fset["eastward_wind"];
  const atlas::Field & vField = fset["northward_wind"];

  const auto modellevels = static_cast<atlas::idx_t>(uField.levels());
  const auto sc = atlas::functionspace::StructuredColumns(uField.functionspace());

  atlas::Field uvgp = sc.createField<double>(atlas::option::name("uv_gp") |
                                             atlas::option::variables(2) |
                                             atlas::option::levels(modellevels));

  const auto uView = atlas::array::make_view<double, 2>(uField);
  const auto vView = atlas::array::make_view<double, 2>(vField);
  auto uvView = atlas::array::make_view<double, 3>(uvgp);

  for (atlas::idx_t jn = 0; jn < uvView.shape()[0]; ++jn) {
    for (atlas::idx_t jl = 0; jl < uvView.shape()[1]; ++jl) {
      uvView(jn, jl, 0) = uView(jn, jl);
      uvView(jn, jl, 1) = vView(jn, jl);
    }
  }

  uvgp.haloExchange();

  return uvgp;
}
// -----------------------------------------------------------------------------

// Create inner variables from outer variables.
// Excludes meridional and zonal winds.
oops::Variables createInnerVars(const oops::Variables & activeVars) {
  oops::Variables innerVars;

  for (const auto & var : activeVars.variables()) {
    if (var != "eastward_wind" && var != "northward_wind") {
      innerVars.push_back(var);
    }
  }

  return innerVars;
}

// -----------------------------------------------------------------------------

void applyNtimesNplus1SpectralScaling(const oops::Variables & innerNames,
                                      const oops::Variables & outerNames,
                                      const atlas::functionspace::Spectral & specFS,
                                      const atlas::idx_t & totalWavenumber,
                                      atlas::FieldSet & fSet,
                                      const bool inverse = false) {
  const auto zonal_wavenumbers = specFS.zonal_wavenumbers();
  const int nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // copy fields that are not associated with innerNames
  atlas::FieldSet fsetTemp;
  for (const atlas::Field & f : fSet) {
    if (!(innerNames.has(f.name()))) {
      fsetTemp.add(f);
    }
  }

  const double earthRadius = atlas::util::Earth::radius();
  const double squaredEarthRadius = earthRadius * earthRadius;
  atlas::FieldSet fsetScaled;
  for (std::size_t var = 0; var < innerNames.variables().size(); ++var) {
    atlas::Field scaledFld = fSet[innerNames[var]];
    auto fldView = atlas::array::make_view<double, 2>(scaledFld);

    int i(0);
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(totalWavenumber); ++n1) {
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (atlas::idx_t jl = 0; jl < fSet[innerNames[var]].levels(); ++jl) {
            if (inverse) {
              if (n1 != 0) {
                fldView(i, jl) /=  n1 * (n1 + 1) / squaredEarthRadius;
              } else {  // Inverse not defined for total wavenumber 0, assume zero mean.
                fldView(i, jl) = 0.0;
              }
            } else {
              fldView(i, jl) *=  n1 * (n1 + 1) / squaredEarthRadius;
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

static SaberOuterBlockMaker<SpectralToGauss> makerSpectralToGauss_("spectral to gauss");

// -----------------------------------------------------------------------------

SpectralToGauss::SpectralToGauss(const oops::GeometryData & outerGeometryData,
                                 const std::vector<size_t> & activeVariableSizes,
                                 const oops::Variables & outerVars,
                                 const eckit::Configuration & covarConf,
                                 const Parameters_ & params,
                                 const atlas::FieldSet & xb,
                                 const atlas::FieldSet & fg,
                                 const util::DateTime & validTimeOfXbFg)
  : SaberOuterBlockBase(params),
    params_(params),
    activeVars_(params.activeVariables.value().get_value_or(outerVars)),
    innerVars_(createInnerVars(activeVars_)),
    outerVars_(outerVars),
    activeVariableSizes_(activeVariableSizes),
    useWindTransform(outerVars_.has("eastward_wind") && outerVars_.has("northward_wind")),
    gaussFunctionSpace_(outerGeometryData.functionSpace()),
    specFunctionSpace_(2 * atlas::GaussianGrid(gaussFunctionSpace_.grid()).N() - 1),
    trans_(gaussFunctionSpace_, specFunctionSpace_),
    innerGeometryData_(atlas::FunctionSpace(specFunctionSpace_),
                       outerGeometryData.fieldSet(),
                       outerGeometryData.levelsAreTopDown(), outerGeometryData.comm())

{
  oops::Log::trace() << classname() << "::SpectralToGauss starting" << std::endl;
  oops::Log::trace() << classname() << "::SpectralToGauss done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyVectorFields(atlas::FieldSet & spectralWindFieldSet,
                                           atlas::FieldSet & outFieldSet) const {
  // 1- Convert stream function and potential vorticity to divergence and vorticity.
  oops::Log::trace() << classname() << "::multiplyVectorFields starting" << std::endl;
  // Scale by n(n+1) / squaredEarthRadius and rename fields to vorticity and divergence
  if (spectralWindFieldSet.has("streamfunction") &&
      spectralWindFieldSet.has("velocity_potential")) {
    ASSERT(innerVars_.has("streamfunction") && innerVars_.has("velocity_potential"));
    const int N = specFunctionSpace_.truncation();
    applyNtimesNplus1SpectralScaling(oops::Variables({"streamfunction", "velocity_potential"}),
                                     oops::Variables({"vorticity", "divergence"}),
                                     specFunctionSpace_, N, spectralWindFieldSet);
  }

  // 2- Convert divergence and vorticity to eastward and northward wind.
  ASSERT(spectralWindFieldSet.has("divergence") && spectralWindFieldSet.has("vorticity"));
  atlas::Field uvgp = allocateGaussUVField(gaussFunctionSpace_,
                                           innerVars_, activeVariableSizes_);
  // Transform to Gaussian grid
  trans_.invtrans_vordiv2wind(spectralWindFieldSet["vorticity"],
                              spectralWindFieldSet["divergence"],
                              uvgp);

  const atlas::FieldSet uvfset = convertUVToFieldSet(uvgp);

  ASSERT(!outFieldSet.has("eastward_wind") && !outFieldSet.has("northward_wind"));
  outFieldSet.add(uvfset["eastward_wind"]);
  outFieldSet.add(uvfset["northward_wind"]);
  oops::Log::trace() << classname() << "::multiplyVectorFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyVectorFieldsAD(atlas::FieldSet & windFieldSet,
                                             atlas::FieldSet & outFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyVectorFieldsAD starting" << std::endl;
  ASSERT(windFieldSet.has("eastward_wind") && windFieldSet.has("northward_wind"));
  ASSERT(!outFieldSet.has("vorticity") && !outFieldSet.has("divergence"));
  ASSERT(!outFieldSet.has("streamfunction") && !outFieldSet.has("velocity_potential"));

  atlas::Field uvgp = convertUVToFieldSetAD(windFieldSet);

  atlas::FieldSet spectralWindFieldSet = allocateSpectralVortDiv(specFunctionSpace_,
                                                                 innerVars_,
                                                                 activeVariableSizes_);

  trans_.invtrans_vordiv2wind_adj(uvgp,
                                  spectralWindFieldSet["vorticity"],
                                  spectralWindFieldSet["divergence"]);


  if (innerVars_.has("streamfunction") && innerVars_.has("velocity_potential")) {
    const int N = specFunctionSpace_.truncation();
    applyNtimesNplus1SpectralScaling(oops::Variables({"vorticity", "divergence"}),
                                     oops::Variables({"streamfunction", "velocity_potential"}),
                                     specFunctionSpace_, N, spectralWindFieldSet);
  }

  outFieldSet.add(spectralWindFieldSet[0]);
  outFieldSet.add(spectralWindFieldSet[1]);
  oops::Log::trace() << classname() << "::multiplyVectorFieldsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyScalarFields(const atlas::FieldSet & specFieldSet,
                                           atlas::FieldSet & outFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyScalarFields starting" << std::endl;
  // Create fieldset on Gaussian grid
  atlas::FieldSet gaussFieldSet;
  for (const auto & fieldname : specFieldSet.field_names()) {
      atlas::Field gaussField =
        gaussFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                   atlas::option::levels(specFieldSet[fieldname].levels()));
      gaussField.haloExchange();
      atlas::array::make_view<double, 2>(gaussField).assign(0.0);
      gaussFieldSet.add(gaussField);
  }

  // Transform to gaussian grid
  trans_.invtrans(specFieldSet, gaussFieldSet);

  for (const auto & fieldname : specFieldSet.field_names()) {
    gaussFieldSet[fieldname].haloExchange();
    ASSERT(!outFieldSet.has(fieldname));
    outFieldSet.add(gaussFieldSet[fieldname]);
  }
  oops::Log::trace() << classname() << "::multiplyScalarFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyScalarFieldsAD(const atlas::FieldSet & gaussFieldSet,
                                             atlas::FieldSet & outFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyScalarFieldsAD starting" << std::endl;
  // Create spectral fieldset
  atlas::FieldSet specFieldSet;
  for (const auto & fieldname : gaussFieldSet.field_names()) {
    atlas::Field specField =
      specFunctionSpace_.createField<double>(atlas::option::name(fieldname) |
                                 atlas::option::levels(gaussFieldSet[fieldname].levels()));
    specFieldSet.add(specField);
  }

  // Transform to spectral space
  trans_.invtrans_adj(gaussFieldSet, specFieldSet);

  for (const auto & fieldname : gaussFieldSet.field_names()) {
    outFieldSet.add(specFieldSet[fieldname]);
  }
  oops::Log::trace() << classname() << "::multiplyScalarFieldsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::invertMultiplyScalarFields(const atlas::FieldSet & gaussFieldSet,
                                                 atlas::FieldSet & outFieldSet) const {
  oops::Log::trace() << classname() << "::invertMultiplyScalarFields starting" << std::endl;

  for (const auto & gaussField : gaussFieldSet) {
    const auto fieldConfig = atlas::option::name(gaussField.name()) |
                             atlas::option::levels(gaussField.levels());
    auto spectralField = specFunctionSpace_.createField<double>(fieldConfig);
    trans_.dirtrans(gaussField, spectralField);
    ASSERT(!outFieldSet.has(spectralField.name()));
    outFieldSet.add(spectralField);
  }

  oops::Log::trace() << classname() << "::invertMultiplyScalarFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::invertMultiplyVectorFields(const atlas::FieldSet & gaussFieldSet,
                                                 atlas::FieldSet & outFieldSet) const {
  oops::Log::trace() << classname() << "::invertMultiplyVectorFields starting" << std::endl;

  ASSERT(gaussFieldSet.has("eastward_wind") && gaussFieldSet.has("northward_wind"));
  ASSERT(!outFieldSet.has("vorticity") && !outFieldSet.has("divergence"));
  ASSERT(!outFieldSet.has("streamfunction") && !outFieldSet.has("velocity_potential"));

  atlas::Field uvgp = convertFieldSetToUV(gaussFieldSet);

  atlas::FieldSet spectralWindFieldSet = allocateSpectralVortDiv(specFunctionSpace_,
                                                                 innerVars_,
                                                                 activeVariableSizes_);

  trans_.dirtrans_wind2vordiv(uvgp,
                              spectralWindFieldSet["vorticity"],
                              spectralWindFieldSet["divergence"]);

  if (innerVars_.has("streamfunction") && innerVars_.has("velocity_potential")) {
    const int N = specFunctionSpace_.truncation();
    const bool inverse = true;
    applyNtimesNplus1SpectralScaling(oops::Variables({"vorticity", "divergence"}),
                                     oops::Variables({"streamfunction", "velocity_potential"}),
                                     specFunctionSpace_, N, spectralWindFieldSet,
                                     inverse);
  }

  outFieldSet.add(spectralWindFieldSet[0]);
  outFieldSet.add(spectralWindFieldSet[1]);

  oops::Log::trace() << classname() << "::invertMultiplyVectorFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet spectralFieldSet = atlas::FieldSet();
  atlas::FieldSet spectralWindFieldSet = atlas::FieldSet();

  // copy "passive variables" and sort out wind variables.
  for (const auto & fieldname : fieldSet.field_names()) {
    if (!activeVars_.has(fieldname)) {  // Passive variables
      newFields.add(fieldSet[fieldname]);
    } else if (useWindTransform && (fieldname == "vorticity" ||
                                    fieldname == "divergence" ||
                                    fieldname == "streamfunction" ||
                                    fieldname == "velocity_potential")) {
        // Active variables to be converted to vector wind
        spectralWindFieldSet.add(fieldSet[fieldname]);
    } else if (fieldname != "eastward_wind" && fieldname != "northward_wind") {
        // Active scalar variables
        spectralFieldSet.add(fieldSet[fieldname]);
    }
  }

  // Convert active wind variables to u/v on Gaussian grid.
  if (useWindTransform) {multiplyVectorFields(spectralWindFieldSet, newFields);}

  // Convert active scalar variables to Gaussian grid.
  if (!spectralFieldSet.empty()) {multiplyScalarFields(spectralFieldSet, newFields);}

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::multiplyAD(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // On input: fieldset on gaussian grid
  atlas::FieldSet newFields = atlas::FieldSet();
  atlas::FieldSet gaussFieldSet = atlas::FieldSet();
  atlas::FieldSet windFieldSet = atlas::FieldSet();

  // copy "passive variables" and sort out wind variables.
  for (const auto & fieldname : fieldSet.field_names()) {
    if (!activeVars_.has(fieldname)) {  // Passive variables
      newFields.add(fieldSet[fieldname]);
    } else if (fieldname == "eastward_wind" || fieldname == "northward_wind") {
        // Active vector wind variables
        windFieldSet.add(fieldSet[fieldname]);
    } else if (!useWindTransform || (fieldname != "vorticity" &&
                                     fieldname != "divergence" &&
                                     fieldname != "streamfunction" &&
                                     fieldname != "velocity_potential")) {
        // Active scalar variables
        gaussFieldSet.add(fieldSet[fieldname]);
    }  // end if active variables
  }  // end for

  // Convert active scalar variables to spectral space.
  if (!gaussFieldSet.empty()) {multiplyScalarFieldsAD(gaussFieldSet, newFields);}

  // Convert active u/v wind variables to spectral space.
  if (useWindTransform) {multiplyVectorFieldsAD(windFieldSet, newFields);}

  fieldSet = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::leftInverseMultiply(atlas::FieldSet & fieldSet) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;
  auto outFieldSet = atlas::FieldSet();
  auto scalarFieldSet = atlas::FieldSet();
  auto windFieldSet = atlas::FieldSet();

  for (const auto & fieldName : fieldSet.field_names()) {
    if (!activeVars_.has(fieldName)) {
        outFieldSet.add(fieldSet[fieldName]);
    } else if (fieldName == "eastward_wind" || fieldName == "northward_wind") {
        windFieldSet.add(fieldSet[fieldName]);
    } else if (!useWindTransform || (fieldName != "vorticity" &&
                                     fieldName != "divergence" &&
                                     fieldName != "streamfunction" &&
                                     fieldName != "velocity_potential")) {
        scalarFieldSet.add(fieldSet[fieldName]);
    }
  }

  if (!scalarFieldSet.empty()) invertMultiplyScalarFields(scalarFieldSet, outFieldSet);

  if (useWindTransform) invertMultiplyVectorFields(windFieldSet, outFieldSet);

  fieldSet = outFieldSet;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SpectralToGauss::print(std::ostream & os) const {
  os << classname();
}

}  // namespace spectralb
}  // namespace saber
