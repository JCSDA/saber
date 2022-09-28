/*
 * (C) Crown Copyright 2020-2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/util/Earth.h"

#include "oops/util/Logger.h"

#include "saber/spectralb/gaussutils.h"

using atlas::grid::detail::partitioner::TransPartitioner;

namespace saber {

std::vector<std::size_t> createActiveVariableSizes(const oops::Variables & activeVars,
                                                   const oops::Variables & inputVars,
                                                   const std::vector<std::size_t> & variableSizes) {
  std::vector<std::size_t> activeVariableSizes(activeVars.size(), 0.0);
  for (std::size_t i = 0; i < activeVars.variables().size(); ++i) {
    activeVariableSizes[i] = variableSizes[inputVars.find(activeVars[i])];
  }
  return activeVariableSizes;
}

void applyNtimesNplus1SpectralScaling(const oops::Variables & inputNames,
                                      const oops::Variables & outputNames,
                                      const atlas::functionspace::Spectral & specFS,
                                      const atlas::idx_t & totalWavenumber,
                                      atlas::FieldSet & fSet) {

  const auto zonal_wavenumbers = specFS.zonal_wavenumbers();
  const int nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // copy fields that are not associated with inputNames
  atlas::FieldSet fsetTemp;
  for (atlas::Field & f : fSet) {
    if (!(inputNames.has(f.name()))) {
      fsetTemp.add(f);
    }
  }

  std::cout << "scaling fsetTemp size" << fsetTemp.size() << std::endl;

  atlas::FieldSet fsetScaled;
  for (std::size_t var = 0; var < inputNames.variables().size(); ++var) {
    atlas::Field scaledFld = fSet[inputNames[var]];
    auto fldView = atlas::array::make_view<double, 2>(scaledFld);

    double a = atlas::util::Earth::radius();  // radius of earth;
    int i(0);
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= static_cast<std::size_t>(2 * totalWavenumber - 1); ++n1) {
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          std::cout << "scaling i n1 m1" << i << " " << n1 << " " << m1 << std::endl;
          for (atlas::idx_t jl = 0; jl < fSet[inputNames[var]].levels(); ++jl) {
            fldView(i, jl) = n1 * (n1 + 1) *  fldView(i, jl) / a;
          }
        }
      }
    }
    scaledFld.rename(outputNames[var]);
    fsetScaled.add(scaledFld);
  }

  fSet = fsetTemp;

  for (auto & f : fsetScaled) {
    fSet.add(f);
  }

  for (auto & f : fSet) {
    std::cout<< "scaling names " << f.name()  << std::endl;
    auto fldView = atlas::array::make_view<double, 2>(f);
    for (atlas::idx_t jnode = 0; jnode < f.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < f.shape(1); ++jlevel) {
        std::cout << f.name() << " " << jnode << " " << jlevel << fldView(jnode, jlevel)
                  <<  " " << fldView(jnode, jlevel) * atlas::util::Earth::radius() << std::endl;
      }
    }


  }

}

atlas::functionspace::Spectral
    createSpectralFunctionSpace(const atlas::StructuredGrid & gaussGrid,
                                const std::vector<std::size_t> & variableSizes) {
  auto N = atlas::GaussianGrid(gaussGrid).N();
  return  atlas::functionspace::Spectral(2*N-1,
       atlas::option::levels(static_cast<atlas::idx_t>(variableSizes[0])));
}

atlas::Field allocateGaussUVField(
    const atlas::FunctionSpace & gaussFS,
    const oops::Variables & activeVariables,
    const std::vector<std::size_t> & activeVariableSizes) {

  std::cout << "start " << activeVariableSizes.size() << std::endl;
  for (std::size_t i : activeVariableSizes) {
    std::cout << i << std::endl;
  }

  for (auto s : activeVariables.variables()) {
    std::cout << s << std::endl;
  }

  std::array<size_t, 2> indx{{0,0}};
  std::array<size_t, 2> levels{{0,0}};
  if (activeVariables.has("vorticity") && activeVariables.has("divergence")) {
    indx[0] = activeVariables.find("vorticity");
    indx[1] = activeVariables.find("divergence");
  } else if (activeVariables.has("streamfunction") && activeVariables.has("velocity_potential")) {
    indx[0] = activeVariables.find("streamfunction");
    indx[1] = activeVariables.find("velocity_potential");
  } else {
    // error trap
    oops::Log::error() << "ERROR - either vorticity and divergence "
                       << "or streamfunction and velocity_potential "
                       << "not present " << std::endl;
    throw std::runtime_error("input fields mis-specified");
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

  std::cout << "before model levels" << std::endl;

  atlas::idx_t modellevels = static_cast<atlas::idx_t>(levels[0]);

  auto sc = atlas::functionspace::StructuredColumns(gaussFS);
  atlas::Field uvgp = sc.createField<double>(atlas::option::name("uv_gp") |
                                             atlas::option::variables(2) |
                                             atlas::option::levels(modellevels));
  return uvgp;
}

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
    throw std::runtime_error("input fields mis-specified");
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



atlas::functionspace::StructuredColumns
createGaussFunctionSpace(const atlas::StructuredGrid & gaussGrid) {
  oops::Log::trace() << "inside createGaussFunctionSpace" << std::endl;
  return atlas::functionspace::StructuredColumns(
    gaussGrid,
    atlas::grid::Partitioner(new TransPartitioner()),
    atlas::option::halo(1));
}

atlas::FieldSet allocateGaussFieldset(
    const atlas::functionspace::StructuredColumns & gaussFunctionSpace,
    const oops::Variables & gaussNames,
    const std::vector<std::size_t> & variableSizes) {

  // create gauss FieldSet (fields are zeroed)
  atlas::FieldSet gaussFieldSet("gauss fieldset with trans partitioner");

  std::size_t i(0);
  for (auto var : gaussNames.variables()) {
    atlas::Field gaussField =
      gaussFunctionSpace.createField<double>(
        atlas::option::name(var) |
        atlas::option::levels(
          static_cast<atlas::idx_t>(variableSizes[i])));

    gaussField.haloExchange();
    gaussFieldSet.add(gaussField);
    ++i;
  }

  return gaussFieldSet;
}

}  // namespace saber
