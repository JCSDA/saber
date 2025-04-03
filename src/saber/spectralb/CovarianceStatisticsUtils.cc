/*
 * (C) Crown Copyright 2017-2025 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "saber/spectralb/CovarianceStatisticsUtils.h"

#include <netcdf.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/Logger.h"

#include "saber/spectralb/spectralb_covstats_interface.h"
#include "saber/spectralb/spectralbParameters.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace saber {
namespace spectralb {
namespace specutils {

using atlas::array::make_view;
using atlas::idx_t;

void copySpectralFieldSet(const atlas::FieldSet & otherFset,
                          atlas::FieldSet & fset) {
  ASSERT(fset.empty());
  for (idx_t ivar = 0; ivar < otherFset.size(); ++ivar) {
    auto spectralMatrix =
      atlas::Field(otherFset.field_names()[ivar],
                   atlas::array::make_datatype<double>(),
                   atlas::array::make_shape(otherFset[ivar].shape(0),
                                            otherFset[ivar].shape(1),
                                            otherFset[ivar].shape(2)));
    auto spectralMatrixView = make_view<double, 3>(spectralMatrix);
    spectralMatrixView.assign(0.0);
    fset.add(spectralMatrix);
  }
  for (atlas::Field & field : fset) {
    auto view = make_view<double, 3>(field);
    auto otherView = make_view<const double, 3>(otherFset[field.name()]);
    for (idx_t jn = 0; jn < field.shape(0); ++jn) {
      for (idx_t jl = 0; jl < field.shape(1); ++jl) {
        for (idx_t jl2 = 0; jl2 < field.shape(2); ++jl2) {
          view(jn, jl, jl2) = otherView(jn, jl, jl2);
        }
      }
    }
    // Copy metadata
    field.metadata() = otherFset[field.name()].metadata();
  }
}


atlas::FieldSet createCorrelUMatrices(const oops::Variables & activeVars,
                                      const atlas::FieldSet & spectralVerticalCovariances,
                                      const atlas::FieldSet & spectralUMatrices,
                                      const atlas::FieldSet & verticalSDs) {
  atlas::FieldSet spectralCorrelUMatrices;

  for (std::size_t i = 0; i < activeVars.size(); ++i) {
    std::string var = activeVars[i].name();

    const auto uMatrixView = make_view<double, 3>(spectralUMatrices[var]);
    const auto verticalSDView = make_view<const double, 1>(verticalSDs[var]);
    const idx_t nSpectralBins = spectralVerticalCovariances[var].shape(0);
    const double sqrtNSpectralBins = std::sqrt(static_cast<double>(nSpectralBins));

    auto correlUMatrix = atlas::Field(var,
                                      atlas::array::make_datatype<double>(),
                                      atlas::array::make_shape(uMatrixView.shape(0),
                                                               uMatrixView.shape(1),
                                                               uMatrixView.shape(2)));
    auto correlUMatrixView = make_view<double, 3>(correlUMatrix);

    for (idx_t bin = 0; bin < uMatrixView.shape(0); ++bin) {
      for (idx_t k1 = 0; k1 < uMatrixView.shape(1); ++k1) {
        for (idx_t k2 = 0; k2 < uMatrixView.shape(2); ++k2) {
          correlUMatrixView(bin, k1, k2) = uMatrixView(bin, k1, k2)
                                           * sqrtNSpectralBins / verticalSDView(k1);
        }
      }
    }

    spectralCorrelUMatrices.add(correlUMatrix);
  }

  return spectralCorrelUMatrices;
}


atlas::FieldSet createSpectralCorrelations(const oops::Variables & activeVars,
                                           const atlas::FieldSet & spectralVerticalCovariances,
                                           const atlas::FieldSet & verticalSDs) {
  atlas::FieldSet spectralCorrelations;

  for (const auto & var : activeVars) {
    const int modelLevels = var.getLevels();
    auto spectralVertCovView =
      make_view<const double, 3>(spectralVerticalCovariances[var.name()]);

    auto verticalSDView =
      make_view<const double, 1>(verticalSDs[var.name()]);

    auto spectralVertCorrel = atlas::Field(var.name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(spectralVertCovView.shape(0),
                               spectralVertCovView.shape(1),
                               spectralVertCovView.shape(2)));

    auto correlationScaling = atlas::Field(var.name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels, modelLevels));

    auto spectralVertCorrelView =
      make_view<double, 3>(spectralVertCorrel);

    auto correlationScalingView =
      make_view<double, 2>(correlationScaling);

    const idx_t nSpectralBins = spectralVerticalCovariances[var.name()].shape(0);

    for (int k1 = 0; k1 < modelLevels; ++k1) {
      for (int k2 = 0; k2 < modelLevels; ++k2) {
        // Assumes the total variance per level is distributed equally across spectral bins.
        correlationScalingView(k1, k2) =
          static_cast<double>(nSpectralBins) /
          (verticalSDView(k1) * verticalSDView(k2));
      }
    }

    for (idx_t bin = 0; bin < spectralVertCorrel.shape(0); ++bin) {
      for (idx_t k1 = 0; k1 < spectralVertCorrel.shape(1); ++k1) {
        for (idx_t k2 = 0; k2 < spectralVertCorrel.shape(2); ++k2) {
          spectralVertCorrelView(bin, k1, k2) = spectralVertCovView(bin, k1, k2) *
            correlationScalingView(k1, k2);
        }
      }
    }

    spectralCorrelations.add(spectralVertCorrel);
  }

  return spectralCorrelations;
}


void createSpectralCovarianceFromUMatrixFile(const std::string & var,
                                             const std::string & netcdfvar,
                                             const spectralbReadParameters & readparams,
                                             atlas::Field & vertcov) {
  auto uMatrix = atlas::Field(var, atlas::array::make_datatype<double>(),
    atlas::array::make_shape(vertcov.shape(0), vertcov.shape(1), vertcov.shape(2)));

  // get covlevels and covbins
  int covbins(0);
  int covlevels(0);

  covSpectralBinsLevels_f90(readparams.toConfiguration(),
                            static_cast<int>(netcdfvar.size()),
                            netcdfvar.c_str(),
                            covbins,
                            covlevels);

  // get umatrix
  const int sizeVec = covlevels * covlevels * covbins;

  std::vector<float> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);

  covSpectralUMatrix_f90(readparams.toConfiguration(),
                         static_cast<int>(netcdfvar.size()),
                         netcdfvar.c_str(),
                         covbins,
                         sizeVec,
                         spectralUMatrix1D[0]);

  const int loff = readparams.levelOffset;  // starting from 0

  auto uMatrixView = make_view<double, 3>(uMatrix);
  std::size_t jn(0);
  for (idx_t bin = 0; bin < static_cast<idx_t>(vertcov.shape(0)); ++bin) {
    for (idx_t k1 = 0; k1 < static_cast<idx_t>(covlevels); ++k1) {
      for (idx_t k2 = 0; k2 < static_cast<idx_t>(covlevels); ++k2, ++jn) {
        if ((k1 >= loff) && (k2 >= loff) &&
            (k1 < loff + vertcov.shape(1)) && (k2 < loff + vertcov.shape(2))) {
          uMatrixView(bin, k1 - loff, k2 - loff) = spectralUMatrix1D[jn];
        }
      }
    }
  }

  // calculate vertical covariances
  auto spectralVertCovView = make_view<double, 3>(vertcov);

  double val;
  for (idx_t bin = 0; bin < spectralVertCovView.shape(0); ++bin) {
    for (idx_t k1 = 0; k1 < spectralVertCovView.shape(1); ++k1) {
      for (idx_t k2 = 0; k2 < spectralVertCovView.shape(2); ++k2) {
        val = 0.0;
        for (idx_t k3 = 0; k3 < uMatrixView.shape(2); ++k3) {
          val += uMatrixView(bin, k1, k3) * uMatrixView(bin, k2, k3);
        }
        // Currently the true vertical covariances and power spectra are
        // multiplied by the number of spectral bins.
        // covbins holds the number of spectral bins that are in the cov file.
        // spectralVertCovView.shape(0) is the number of spectral bins that we
        // use, that is consistent with the Gaussian resolution employed.
        // So we rescale the covariances to be consistent with this reduced
        // number of covariance bins (and this contract).
        // TODO(MW) - Investigate whether we can remove this scaling throughout
        //            the code.
        spectralVertCovView(bin, k1, k2) = val *
            static_cast<double>(spectralVertCovView.shape(0)) /
            static_cast<double>(covbins);
      }
    }
  }
}


void createUMatrixFromSpectralCovarianceFile(const std::string & var,
                                             const spectralbReadParameters & readparams,
                                             atlas::Field & umatrix) {
  auto spectralVertCov = atlas::Field(var, atlas::array::make_datatype<double>(),
    atlas::array::make_shape(umatrix.shape(0), umatrix.shape(1), umatrix.shape(2)));

  readSpectralCovarianceFromFile(var, var, readparams, spectralVertCov);

  const int sizeVec = umatrix.shape(1) * umatrix.shape(2);
  std::vector<double> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);
  auto spectralVertCovView = make_view<double, 3>(spectralVertCov);
  auto umatrixView = make_view<double, 3>(umatrix);
  umatrixView.assign(0.0);

  for (idx_t bin = 0; bin < static_cast<idx_t>(umatrix.shape(0)); ++bin) {
    std::size_t jn(0);
    for (idx_t k1 = 0; k1 < static_cast<idx_t>(umatrix.shape(1)); ++k1) {
      for (idx_t k2 = 0; k2 < static_cast<idx_t>(umatrix.shape(2)); ++k2, ++jn) {
        spectralUMatrix1D[jn] = spectralVertCovView(bin, k1, k2);
      }
    }

    calculatingSqrtB_f90(umatrix.shape(1), spectralUMatrix1D[0]);

    jn = 0;
    for (idx_t k1 = 0; k1 < static_cast<idx_t>(umatrix.shape(1)); ++k1) {
      for (idx_t k2 = 0; k2 < static_cast<idx_t>(umatrix.shape(2)); ++k2, ++jn) {
        // index order here switched to k2, k1 because of square root being calculated in
        // Fortran
        umatrixView(bin, k2, k1) =  k1 > k2 ? 0.0  : spectralUMatrix1D[jn];
      }
    }
  }
}


atlas::FieldSet createSpectralCovariances(const oops::Variables & activeVars,
                                          const std::vector<std::size_t> & nSpectralBinsFull,
                                          const std::size_t nSpectralBins,
                                          const atlas::FieldSet & spectralUMatrices)
{
  atlas::FieldSet spectralVerticalCovariances;

  for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
    const std::string var = activeVars[ivar].name();
    const int modelLevels = activeVars[ivar].getLevels();
    const std::size_t nSpectralBinsUse = std::min(nSpectralBins, nSpectralBinsFull[ivar]);

    auto spectralVertCov = atlas::Field(var, atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nSpectralBinsUse, modelLevels, modelLevels));

    auto spectralVertCovView =
      make_view<double, 3>(spectralVertCov);
    auto uMatrixView = make_view<double, 3>(spectralUMatrices[var]);

    double val;
    for (idx_t bin = 0; bin < spectralVertCovView.shape(0); ++bin) {
      for (idx_t k1 = 0; k1 < spectralVertCovView.shape(1); ++k1) {
        for (idx_t k2 = 0; k2 < spectralVertCovView.shape(2); ++k2) {
          val = 0.0;
          for (idx_t k3 = 0; k3 < uMatrixView.shape(2); ++k3) {
            val += uMatrixView(bin, k1, k3) * uMatrixView(bin, k2, k3);
          }
          // There is a loss of variance here when nSpectralBins is less than nSpectralBinsFull,
          // as we only keep nSpectralBins out of
          // nSpectralBinsFull[ivar].A crude renormalization is applied, assuming
          // the variance is equally distributed across bins:
          spectralVertCovView(bin, k1, k2) = val * nSpectralBins / (nSpectralBinsFull[ivar]);
        }
      }
    }

    spectralVerticalCovariances.add(spectralVertCov);
  }

  return spectralVerticalCovariances;
}

atlas::FieldSet createUMatrices(
    const oops::Variables & activeVars,
    const std::vector<std::size_t> & nSpectralBinsFull,
    const spectralbReadParameters & params) {
  const auto & umatrixNetCDFParams = params.umatrixNetCDFNames.value();
  atlas::FieldSet spectralUMatrices;

  if (umatrixNetCDFParams == boost::none) {
    for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
      const std::string var = activeVars[ivar].name();
      const int modelLevels = activeVars[ivar].getLevels();
      auto uMatrix = atlas::Field(var, atlas::array::make_datatype<double>(),
        atlas::array::make_shape(nSpectralBinsFull[ivar], modelLevels, modelLevels));
      createUMatrixFromSpectralCovarianceFile(var,
                                              params,
                                              uMatrix);
      spectralUMatrices.add(uMatrix);
    }
  } else {
    oops::Variables netCDFVars(umatrixNetCDFParams.value());
    for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
      std::string var = activeVars[ivar].name();
      const int modelLevels = activeVars[ivar].getLevels();
      std::string netCDFVar = netCDFVars[ivar].name();

      auto uMatrix = atlas::Field(var, atlas::array::make_datatype<double>(),
        atlas::array::make_shape(nSpectralBinsFull[ivar], modelLevels, modelLevels));

      // vector size
      const int sizeVec = modelLevels * modelLevels * static_cast<int>(nSpectralBinsFull[ivar]);

      std::vector<float> spectralUMatrix1D(static_cast<std::size_t>(sizeVec), 0.0);

      covSpectralUMatrix_f90(params.toConfiguration(),
                             static_cast<int>(netCDFVar.size()),
                             netCDFVar.c_str(),
                             static_cast<int>(nSpectralBinsFull[ivar]),
                             sizeVec,
                             spectralUMatrix1D[0]);


      auto uMatrixView = make_view<double, 3>(uMatrix);
      std::size_t jn(0);
      for (idx_t bin = 0; bin < static_cast<idx_t>(nSpectralBinsFull[ivar]); ++bin) {
        for (idx_t k1 = 0; k1 < static_cast<idx_t>(modelLevels); ++k1) {
          for (idx_t k2 = 0; k2 < static_cast<idx_t>(modelLevels); ++k2, ++jn) {
            uMatrixView(bin, k1, k2) = spectralUMatrix1D[jn];
          }
        }
      }
      spectralUMatrices.add(uMatrix);
    }
  }

  return spectralUMatrices;
}


atlas::FieldSet createVerticalSD(const oops::Variables & activeVars,
                                 const atlas::FieldSet & spectralVerticalCovariances) {
  atlas::FieldSet verticalSDs;

  for (const auto & var : activeVars) {
    const int modelLevels = var.getLevels();
    auto spectralVertCovView =
      make_view<const double, 3>(spectralVerticalCovariances[var.name()]);

    auto verticalSD = atlas::Field(var.name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(modelLevels));

    auto verticalSDView = make_view<double, 1>(verticalSD);

    for (int k = 0; k < modelLevels; ++k) {
      verticalSDView(k) = 0.0;
      for (idx_t bin = 0; bin < spectralVertCovView.shape(0); ++bin) {
        verticalSDView(k) += spectralVertCovView(bin, k, k);
      }
      verticalSDView(k) = std::sqrt(verticalSDView(k));
    }
    verticalSDs.add(verticalSD);
  }

  return verticalSDs;
}

std::vector<std::size_t> getNSpectralBinsFull(const spectralbReadParameters & params,
                                              const oops::Variables & activeVars) {
  const auto & umatrixNetCDFParams = params.umatrixNetCDFNames.value();

  std::vector<std::size_t> nSpectralBinsFull(activeVars.size());

  if (umatrixNetCDFParams == boost::none) {
    // Setup
    const std::string filepath = params.covarianceFile;
    std::vector<std::string> dimNames;
    std::vector<idx_t> dimSizes;
    std::vector<std::vector<std::string>> dimNamesForEveryVar;
    oops::Variables vars;
    std::vector<int> netcdfGeneralIDs;
    eckit::LocalConfiguration netcdfMetaData;
    std::vector<int> netcdfDimIDs;
    std::vector<int> netcdfVarIDs;
    std::vector<std::vector<int>> netcdfDimVarIDs;

    util::atlasArrayInquire(filepath,
                            dimNames,
                            dimSizes,
                            vars,
                            dimNamesForEveryVar,
                            netcdfMetaData,
                            netcdfGeneralIDs,
                            netcdfDimIDs,
                            netcdfVarIDs,
                            netcdfDimVarIDs);

    // TODO(Marek) - extend so that multiple horizontal wavenumbers possible
    std::string hwaveno("binning_index");
    auto ind = std::find(dimNames.begin(),  dimNames.end(), hwaveno);
    if (ind != dimNames.end()) {
      for (std::size_t ivar = 0; ivar < activeVars.size(); ++ivar) {
        nSpectralBinsFull[ivar] = dimSizes[std::distance(dimNames.begin(), ind)];
      }
    }
  } else {
    const oops::Variables vars(umatrixNetCDFParams.value());

    for (std::size_t ivar = 0; ivar < vars.size(); ++ivar) {
      const std::string var = vars[ivar].name();
      int nbins(0);

      // get the number of spectral bins from the cov file
      covSpectralBins_f90(params.toConfiguration(),
                          static_cast<int>(var.size()),
                          var.c_str(),
                          nbins);

      nSpectralBinsFull[ivar] = static_cast<std::size_t>(nbins);
    }
  }

  return nSpectralBinsFull;
}


void readSpectralCovarianceFromFile(const std::string & varname1,
                                    const std::string & varname2,
                                    const spectralbReadParameters & readparams,
                                    atlas::Field & spectralVertCov) {
  std::string ncfilepath = readparams.covarianceFile;
  std::vector<std::string> dimNames;
  std::vector<idx_t> dimSizes;
  oops::Variables vars;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  std::vector<int> netcdfGeneralIDs;
  eckit::LocalConfiguration netCDFConf;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  // read from header file on root PE.
  const std::size_t root = 0;
  std::string filevarname("");
  if (eckit::mpi::comm().rank() == root) {
     util::atlasArrayInquire(ncfilepath,
                             dimNames,
                             dimSizes,
                             vars,
                             dimNamesForEveryVar,
                             netCDFConf,
                             netcdfGeneralIDs,
                             netcdfDimIDs,
                             netcdfVarIDs,
                             netcdfDimVarIDs);
    for (const oops::Variable & var : vars) {
      const std::string var1 =
        util::getAttributeValue<std::string>(netCDFConf, var.name(),
                                             "variable_name_1");
      const std::string var2 =
        util::getAttributeValue<std::string>(netCDFConf, var.name(),
                                             "variable_name_2");
      if ((var1 == varname1) && (var2 == varname2)) {
        filevarname = var.name();
      }
    }
  }

  auto specvertview = make_view<double, 3>(spectralVertCov);
  specvertview.assign(0.0);

  if (eckit::mpi::comm().rank() == root) {
    const std::size_t i = vars.find(filevarname);

    std::vector<idx_t> dimSizesForVar;
    for (const auto & dimName : dimNamesForEveryVar[i]) {
      auto it2 = std::find(dimNames.begin(), dimNames.end(), dimName);
      ASSERT(it2 != dimNames.end());
      std::size_t i2 = std::distance(dimNames.begin(), it2);
      dimSizesForVar.push_back(dimSizes.at(i2));
    }

    if ( (dimSizesForVar[0] != specvertview.shape(0)) ||
         (dimSizesForVar[1] != specvertview.shape(1)) ||
         (dimSizesForVar[2] != specvertview.shape(2)) ) {
      oops::Log::error() <<
        "Dimensions of spectral vertical covariances in file is inconsistent with" <<
        " yaml setup for variable " <<  filevarname << std::endl;
      throw eckit::UserError(
        "Inconsistency in dimension sizes of spectral vertical covariance"
        " between file and assumed array shape", Here());
    }

    util::atlasArrayReadData(netcdfGeneralIDs,
                             netcdfVarIDs[i],
                             specvertview);
  }
  util::scatter<double>(eckit::mpi::comm(), root, specvertview);

  if (eckit::mpi::comm().rank() == root) {
    const int retval = nc_close(netcdfGeneralIDs[0]);
    if (retval) ERR(retval);
  }
}


void spectralVerticalConvolution(const oops::Variables & activeVars,
                                 const atlas::functionspace::Spectral & specFunctionSpace,
                                 const atlas::FieldSet & spectralVerticalStats,
                                 atlas::FieldSet & fieldSet) {
  const size_t N = specFunctionSpace.truncation();

  const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars) {
    idx_t i = 0;
    idx_t levels(fieldSet[var.name()].shape(1));
    auto spfView = make_view<double, 2>(fieldSet[var.name()]);
    auto vertCovView = make_view<const double, 3>(spectralVerticalStats[var.name()]);
    const size_t Ncov = vertCovView.shape()[0] - 1;

    std::vector<double> col(levels), col2(levels);

    if (N > Ncov) {
      // For each total wavenumber n1, perform a 1D convolution with vertical covariances.
      for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
        const idx_t m1 = zonal_wavenumbers(jm);
        for (std::size_t n1 = m1; n1 <= N; ++n1) {
          // Note that img stands for imaginary component and are the
          // odd indices in the first index of the spectral fields.
          for (std::size_t img = 0; img < 2; ++img) {
            // Pre-fill vertical column to be convolved.
            if (n1 <= Ncov) {
              for (idx_t jl = 0; jl < levels; ++jl) {
                col[jl] = spfView(i, jl);
              }
              // The 2*n1+1 factor is there to equally distribute the covariance across
              // the spectral coefficients associated to this total wavenumber.
              const double norm = static_cast<double>((2 * n1 + 1) * vertCovView.shape(0));
              for (idx_t r = 0; r < levels; ++r) {
                col2[r] = 0;
                for (idx_t c = 0; c < levels; ++c) {
                  col2[r] += vertCovView(n1, r, c) * col[c] / norm;
                }
              }
              for  (idx_t jl = 0; jl < levels; ++jl) {
                spfView(i, jl) = col2[jl];
              }
            } else {
              // we are outside specified statistics so we truncate.
              // we may want to do something different later here
              for (idx_t jl = 0; jl < levels; ++jl) {
                spfView(i, jl) = 0.0;
              }
            }
            ++i;
          }
        }
      }
    } else {
      // For each total wavenumber n1, perform a 1D convolution with vertical covariances.
      for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
        const idx_t m1 = zonal_wavenumbers(jm);
        for (std::size_t n1 = m1; n1 <= N; ++n1) {
          // Note that img stands for imaginary component and are the
          // odd indices in the first index of the spectral fields.
          for (std::size_t img = 0; img < 2; ++img) {
            // Pre-fill vertical column to be convolved.
            for (idx_t jl = 0; jl < levels; ++jl) {
              col[jl] = spfView(i, jl);
            }
            // The 2*n1+1 factor is there to equally distribute the covariance across
            // the spectral coefficients associated to this total wavenumber.
            const double norm = static_cast<double>((2 * n1 + 1) * vertCovView.shape(0));
            for (idx_t r = 0; r < levels; ++r) {
              col2[r] = 0;
              for (idx_t c = 0; c < levels; ++c) {
                col2[r] += vertCovView(n1, r, c) * col[c] / norm;
              }
            }
            for  (idx_t jl = 0; jl < levels; ++jl) {
              spfView(i, jl) = col2[jl];
            }
            ++i;
          }
        }
      }
    }
  }
}


void spectralVerticalConvolutionSqrt(const oops::Variables & activeVars,
                                     const atlas::functionspace::Spectral & specFunctionSpace,
                                     const atlas::FieldSet & spectralVerticalStatsSqrt,
                                     atlas::FieldSet & fieldSet) {
  const size_t N = specFunctionSpace.truncation();

  const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars) {
    idx_t levels(fieldSet[var.name()].shape(1));
    auto UMatrixView = make_view<const double, 3>(spectralVerticalStatsSqrt[var.name()]);
    auto spfView = make_view<double, 2>(fieldSet[var.name()]);
    const int nSpectralBinsFull = spectralVerticalStatsSqrt[var.name()].shape(0);

    idx_t i = 0;
    std::vector<double> col(levels), col2(levels);
    for (idx_t jm1 = 0; jm1 < nb_zonal_wavenumbers; ++jm1) {
      const idx_t m1 = zonal_wavenumbers(jm1);
      for (std::size_t n1 = m1; n1 <= N; ++n1) {
        // note that img stands for imaginary component are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[jl] = spfView(i, jl);
          }
          for (idx_t r = 0; r < levels; ++r) {
            col2[r] = 0;
            for (idx_t c = 0; c < levels; ++c) {
              col2[r] += UMatrixView(n1, r, c) * col[c];
            }
          }
          const double norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[jl] / norm;
          }
        }
      }
    }
  }
}


void spectralVerticalConvolutionSqrtAD(const oops::Variables & activeVars,
                                       const atlas::functionspace::Spectral & specFunctionSpace,
                                       const atlas::FieldSet & spectralVerticalStatsSqrt,
                                       atlas::FieldSet & fieldSet) {
  const size_t N = specFunctionSpace.truncation();

  const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  // Only update the fields that were specified in the active variables
  for (const auto & var : activeVars) {
    idx_t levels(fieldSet[var.name()].shape(1));
    auto UMatrixView = make_view<const double, 3>(spectralVerticalStatsSqrt[var.name()]);
    auto spfView = make_view<double, 2>(fieldSet[var.name()]);
    const int nSpectralBinsFull = spectralVerticalStatsSqrt[var.name()].shape(0);

    idx_t i = 0;
    std::vector<double> col(levels), col2(levels);
    for (idx_t jm1 = 0; jm1 < nb_zonal_wavenumbers; ++jm1) {
      const idx_t m1 = zonal_wavenumbers(jm1);
      for (std::size_t n1 = m1; n1 <= N; ++n1) {
        // note that img stands for imaginary component are the
        // odd indices in the first index of the spectral fields.
        for (std::size_t img = 0; img < 2; ++img, ++i) {
          for (idx_t jl = 0; jl < levels; ++jl) {
            col[jl] = spfView(i, jl);
          }
          for (idx_t r = 0; r < levels; ++r) {
            col2[r] = 0;
            for (idx_t c = 0; c < levels; ++c) {
              col2[r] += UMatrixView(n1, c, r) * col[c];
            }
          }
          const double norm = std::sqrt(static_cast<double>((2 * n1 + 1) * nSpectralBinsFull));
          for  (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) = col2[jl] / norm;
          }
        }
      }
    }
  }
}


/// Same as spectralVerticalConvolution, without the vertical convolution
void spectralHorizontalFilter(const oops::Variables & activeVars,
                              const atlas::functionspace::Spectral & specFunctionSpace,
                              const atlas::FieldSet & spectralVerticalStats,
                              atlas::FieldSet & fieldSet) {
  const size_t N = specFunctionSpace.truncation();
  const auto zonal_wavenumbers = specFunctionSpace.zonal_wavenumbers();
  const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

  for (const auto & var : activeVars) {
    const idx_t levels(fieldSet[var.name()].shape(1));
    auto vertCovView = make_view<const double, 3>(spectralVerticalStats[var.name()]);
    auto spfView = make_view<double, 2>(fieldSet[var.name()]);

    idx_t i = 0;
    for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const idx_t m1 = zonal_wavenumbers(jm);
      for (std::size_t n1 = m1; n1 <= N; ++n1) {
        for (std::size_t img = 0; img < 2; ++img) {
          const double norm = static_cast<double>((2 * n1 + 1) * vertCovView.shape(0));
          for (idx_t jl = 0; jl < levels; ++jl) {
            spfView(i, jl) *= vertCovView(n1, jl, jl) / norm;
          }
          ++i;
        }
      }
    }
  }
}



std::size_t updateSpectralVerticalCovariances(
    const oops::FieldSets & ensFieldSet,
    const std::size_t priorSampleSize,
    atlas::FieldSet & spectralVerticalCovariances) {

  // assuming that the fields in each fieldset are consistent across the
  // ensemble (ie the same size and the same type of functionspace)

  // we expect with priorSampleSize = 0 for the vertical covariances to be zero.

  // we expect that (for now) that total wavenumbers start from 0
  // TO DO:(Marek) - maybe add metadata to field to include spectral offset?
  // TO DO:(Marek) Normalisation here 1/N (if coming from randomisation).  1/(n-1) if not
  //       (for now set to randomisation version)

  if (priorSampleSize > 0) {
    // assuming read done before and multiplying by number of prior samples
    for (atlas::Field & vertcov : spectralVerticalCovariances) {
      auto vertCovView = make_view<double, 3>(vertcov);
      for (idx_t i = 0; i < vertcov.shape()[0]; ++i) {
        for (idx_t j = 0; j < vertcov.shape()[1]; ++j) {
          for (idx_t k = 0; k < vertcov.shape()[2]; ++k) {
            vertCovView(i, j, k) *= static_cast<double>(priorSampleSize);
          }
        }
      }
    }
  }


  // no read - but instead accumulate matrix only
  for (atlas::Field & vertcov : spectralVerticalCovariances) {
    auto vertCovView = make_view<double, 3>(vertcov);
    int levels = vertcov.shape()[2];
    std::string name(vertcov.name());
    const auto zonal_wavenumbers =
      atlas::functionspace::Spectral(
        ensFieldSet[0].fieldSet()[name].functionspace()).zonal_wavenumbers();
    const size_t N =
      atlas::functionspace::Spectral(ensFieldSet[0].fieldSet()[name].functionspace()).truncation();
    const idx_t nb_zonal_wavenumbers = zonal_wavenumbers.size();

    for (size_t jj = 0; jj < ensFieldSet.size(); ++jj) {
      const auto & fs = ensFieldSet[jj];
      auto spfView = make_view<const double, 2>(fs[name]);

      idx_t i = 0;  // i is the spectral coefficient index
                           // for a given level.
      // For each total wavenumber n1, perform a 1D convolution with vertical covariances.
      // We are assuming that the mean has been removed from the perturbations.
      for (idx_t jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
        const idx_t m1 = zonal_wavenumbers(jm);
        for (std::size_t n1 = m1; n1 <= N; ++n1) {
          // Note 1: that img stands for imaginary component and are the
          // odd indices in the first index of the spectral fields.
          // In the MetOffice system the scaling was 2 * (n_1+1) * n_bins
          // Here the scaling is different ... we scale by n_bins and
          // double count the entries for m/=0 to take account of the
          // implicit complex conjugate spectral coefficients.

          const double scaling = static_cast<double>(vertcov.shape()[0]);
          for (std::size_t img = 0; img < 2; ++img) {
            for (idx_t r = 0; r < levels; ++r) {
              for (idx_t c = 0; c < levels; ++c) {
                if (m1 == 0) {
                  vertCovView(n1, r, c) += spfView(i, r) * spfView(i, c) * scaling;
                } else {
                  vertCovView(n1, r, c) += 2.0 * spfView(i, r) * spfView(i, c) * scaling;
                }
              }
            }
            ++i;
          }
        }
      }
    }
  }

  const std::size_t updatedSampleSize = priorSampleSize + ensFieldSet.size();

  const double recipPriorSampleSize = 1.0/static_cast<double>(updatedSampleSize);
  for (atlas::Field & vertcov : spectralVerticalCovariances) {
    auto vertCovView = make_view<double, 3>(vertcov);
    for (idx_t i = 0; i < vertcov.shape()[0]; ++i) {
      for (idx_t j = 0; j < vertcov.shape()[1]; ++j) {
        for (idx_t k = 0; k < vertcov.shape()[2]; ++k) {
          vertCovView(i, j, k) *= recipPriorSampleSize;
        }
      }
    }
  }
  return  updatedSampleSize;
}

}  // namespace specutils
}  // namespace spectralb
}  // namespace saber
