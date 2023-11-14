/*
 * (C) Crown Copyright 2017-2022 Met Office
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

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"

#include "saber/spectralb/spectralbParameters.h"

namespace saber {
namespace spectralb {

std::vector<std::size_t> getNSpectralBinsFull(const spectralbParameters &);

atlas::FieldSet createUMatrices(const oops::Variables &,
                                const int,
                                const std::vector<std::size_t> &,
                                const spectralbParameters &);

atlas::FieldSet createSpectralCovariances(const oops::Variables &,
                                          const int,
                                          const std::vector<std::size_t> &,
                                          const atlas::FieldSet &,
                                          const spectralbParameters &);

atlas::FieldSet createVerticalSD(const oops::Variables &,
                                 const atlas::FieldSet &);

atlas::FieldSet createCorrelUMatrices(const oops::Variables &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &);

atlas::FieldSet createSpectralCorrelations(const oops::Variables &,
                                           const atlas::FieldSet &,
                                           const atlas::FieldSet &);

void createSpectralCovarianceFromUMatrixFile(const std::string &,
                                             const std::string &,
                                             const spectralbReadVertCovParameters &,
                                             atlas::Field &);

void readSpectralCovarianceFromFile(const std::string &,
                                    const spectralbReadVertCovParameters &,
                                    atlas::Field &);

void updateSpectralVerticalCovariances(
    const std::vector<atlas::FieldSet> & ensFieldSet,
    int & priorSampleSize,
    atlas::FieldSet & spectralVerticalCovariances);

void copySpectralFieldSet(const atlas::FieldSet &,
                          atlas::FieldSet &);

void gatherSumSpectralFieldSet(const eckit::mpi::Comm &,
                               const std::size_t & root,
                               atlas::FieldSet & spectralVerticalCovariances);
}  // namespace spectralb
}  // namespace saber
