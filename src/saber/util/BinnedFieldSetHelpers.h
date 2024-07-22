/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

namespace util {

atlas::FieldSet gatherSumFieldSet(const eckit::mpi::Comm &,
                                  const int root,
                                  const atlas::FieldSet &);

atlas::FieldSet gatherSumFieldSet(const eckit::mpi::Comm &,
                                  const int root,
                                  const std::size_t totalbins,
                                  const atlas::Field & glBinIdx,
                                  const atlas::FieldSet &);

/// \brief Binned field gathering over bins and levels on root PE
atlas::Field gatherSumBinsLevels(
  const eckit::mpi::Comm & comm, const int root,
  const std::size_t & nbins, const std::size_t & levels,
  const atlas::Field & glBinIdx, const atlas::Field & fld);

/// \brief Binned field gathering over bins, levels1 and levels2
///        on the root PE
atlas::Field gatherSumBinsLevelsLevels(
  const eckit::mpi::Comm & comm, const int root,
  const std::size_t & nbins, const std::size_t & levels1,
  const std::size_t & levels2, const atlas::Field & glBinIdx,
  const atlas::Field & fld);

}  // namespace util
