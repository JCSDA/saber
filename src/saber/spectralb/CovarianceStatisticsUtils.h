/*
 * (C) Crown Copyright 2017-2023 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSets.h"
#include "oops/base/Variables.h"

#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {
namespace specutils {

void copySpectralFieldSet(const atlas::FieldSet &,
                          atlas::FieldSet &);

atlas::FieldSet createCorrelUMatrices(const oops::Variables &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &);

atlas::FieldSet createSpectralCorrelations(const oops::Variables &,
                                           const atlas::FieldSet &,
                                           const atlas::FieldSet &);

void createSpectralCovarianceFromUMatrixFile(const std::string &,
                                             const std::string &,
                                             const spectralbReadParameters &,
                                             atlas::Field &);

atlas::FieldSet createSpectralCovariances(const oops::Variables &,
                                          const std::vector<std::size_t> &,
                                          const std::size_t,
                                          const atlas::FieldSet &);

atlas::FieldSet createUMatrices(const oops::Variables &,
                                const std::vector<std::size_t> &,
                                const spectralbReadParameters &);

atlas::FieldSet createVerticalSD(const oops::Variables &,
                                 const atlas::FieldSet &);

void gatherSumSpectralFieldSet(const eckit::mpi::Comm &,
                               const std::size_t root,
                               atlas::FieldSet & spectralVerticalCovariances);

std::vector<std::size_t> getNSpectralBinsFull(const spectralbReadParameters &);

void readSpectralCovarianceFromFile(const std::string &,
                                    const spectralbReadParameters &,
                                    atlas::Field &);

void spectralVerticalConvolution(const oops::Variables &,
                                 const atlas::functionspace::Spectral &,
                                 const atlas::FieldSet &,
                                 atlas::FieldSet &);

void spectralVerticalConvolutionSqrt(const oops::Variables &,
                                     const atlas::functionspace::Spectral &,
                                     const atlas::FieldSet &,
                                     atlas::FieldSet &);

void spectralVerticalConvolutionSqrtAD(const oops::Variables &,
                                       const atlas::functionspace::Spectral &,
                                       const atlas::FieldSet &,
                                       atlas::FieldSet &);

void updateSpectralVerticalCovariances(
  const oops::FieldSets & ensFieldSet,
  int & priorSampleSize,
  atlas::FieldSet & spectralVerticalCovariances);

}  // namespace specutils
}  // namespace spectralb
}  // namespace saber
