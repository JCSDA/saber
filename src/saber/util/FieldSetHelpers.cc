/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <tuple>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Earth.h"

#include "eckit/geometry/Point2.h"
#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "saber/util/FieldSetHelpers.h"

namespace util {
// -----------------------------------------------------------------------------

std::tuple<std::vector<std::vector<double>>,
           std::vector<std::vector<double>>,
           std::vector<size_t>,
           std::vector<size_t>>
sortBySeparationDistance(const eckit::mpi::Comm & comm,
                         const atlas::FunctionSpace & fspace,
                         const atlas::FieldSet & dataFset,
                         const atlas::FieldSet & diracFset,
                         const double maxLength,
                         const bool removeDuplicates) {
  oops::Log::trace() << "sortBySeparationDistance starting" << std::endl;

  // Pull out values of lons, lats, levs where diagFset is unity.
  // Values are gathered on root MPI task w.r.t. geometry communicator.
  const eckit::mpi::Comm & timeComm = oops::mpi::myself();
  auto[lons, lats, levs, subWindows, vals, fieldIndex]
          = util::extractUnityPoints(timeComm, comm, fspace, dataFset, diracFset);

  // Silence "unused variables" compiler warning by casting to void (does nothing)
  (void)subWindows;
  (void)vals;

  // Issue warning if there are several Dirac points in a given level.
  const auto & names = diracFset.field_names();
  for (size_t iDirac = 0; iDirac < lons.size(); ++iDirac) {
    for (size_t jDirac = 0; jDirac < iDirac; ++jDirac) {
      if ((levs[iDirac] == levs[jDirac]) && (fieldIndex[iDirac] == fieldIndex[jDirac])) {
        oops::Log::warning() << "   WARNING!   I found two points at level "
                             << levs[iDirac] << " for field "
                             << names[fieldIndex[iDirac]] << std::endl;
        oops::Log::warning() << "              In most situations, this will give weird results,"
                             << " as the presence of the other Dirac points may make"
                             << " the field non-invariant by rotation around the Dirac points."
                             << std::endl;
      }
    }
  }

  // Sort by variable name and level
  std::vector<size_t> diracIndexes(lons.size());
  size_t counter = 0;
  for (const auto & name : names) {
    const size_t levels = diracFset[name].levels();
    for (size_t lev = 0; lev < levels; lev++) {
      for (size_t iDirac = 0; iDirac < lons.size(); ++iDirac) {
        if ((levs[iDirac] == lev) && (names[fieldIndex[iDirac]] == name)) {
          diracIndexes[counter++] = iDirac;
        }
      }
    }
  }
  ASSERT(counter == lons.size());

  auto applyIndexingDouble = [](const std::vector<size_t> & indexes,
                                std::vector<double> & vec) {
    std::vector<double> sorted(indexes.size(), 0.0);
    std::transform(indexes.cbegin(), indexes.cend(), sorted.begin(),
                   [&](const size_t index){return vec[index];});
    vec = sorted;
  };
  auto applyIndexingSizet = [](const std::vector<size_t> & indexes,
                               std::vector<size_t> & vec) {
    std::vector<size_t> sorted(indexes.size(), 0.0);
    std::transform(indexes.cbegin(), indexes.cend(), sorted.begin(),
                   [&](const size_t index){return vec[index];});
    vec = sorted;
  };
  applyIndexingDouble(diracIndexes, lons);
  applyIndexingDouble(diracIndexes, lats);
  applyIndexingSizet(diracIndexes, levs);
  applyIndexingSizet(diracIndexes, fieldIndex);

  // Broadcast lon, lat, lev & fieldIndex (which are small vectors),
  // to make them accessible from all local atlas Field.
  size_t nDiracs(lons.size());
  comm.broadcast(nDiracs, 0);
  lons.resize(nDiracs);
  lats.resize(nDiracs);
  levs.resize(nDiracs);
  fieldIndex.resize(nDiracs);
  comm.broadcast(lons, 0);
  comm.broadcast(lats, 0);
  comm.broadcast(levs, 0);
  comm.broadcast(fieldIndex, 0);

  // Loop over fields to save (distance, value) pairs
  const auto & ghost = fspace.ghost();
  const auto ghview  = atlas::array::make_view<int, 1>(ghost);
  const auto lonlatView = atlas::array::make_view<double, 2>(fspace.lonlat());
  std::vector<std::vector<double>> locDistances(lons.size());
  std::vector<std::vector<double>> locValues(lons.size());
  for (const auto & name : names) {
    ASSERT(dataFset[name].rank() == 2);
    const auto view = atlas::array::make_view<double, 2>(dataFset[name]);
    for (size_t iDirac = 0; iDirac < nDiracs; iDirac++) {
      if (name == names[fieldIndex[iDirac]]) {
        const auto refPoint = atlas::Point2(lons[iDirac], lats[iDirac]);
        for (atlas::idx_t jnode = 0; jnode < view.shape(0); jnode++) {
          if (!ghview(jnode)) {
            const auto point = atlas::Point2(lonlatView(jnode, 0), lonlatView(jnode, 1));
            const double dist = atlas::util::Earth().distance(refPoint, point);
            if (dist <= maxLength) {
              const double value = view(jnode, levs[iDirac]);
              locDistances[iDirac].push_back(dist);
              locValues[iDirac].push_back(value);
            }
          }
        }
      }
    }
  }

  // Gather distances and values onto root MPI task for sorting
  std::vector<std::vector<double>> distances(nDiracs);
  std::vector<std::vector<double>> values(nDiracs);
  for (size_t iDirac = 0; iDirac < nDiracs; iDirac++) {
    // Gather sizes
    const int dSize = locDistances[iDirac].size();
    std::vector<int> dSizes(comm.size());
    comm.gather(dSize, dSizes, 0);

    // Gather data
    std::vector<int> dDispls;
    std::vector<int> dRecvcounts;
    if (comm.rank() == 0) {
      dDispls.resize(comm.size());
      dRecvcounts.resize(comm.size());
      for (size_t i = 0; i < comm.size(); ++i) {
        dRecvcounts[i] = dSizes[i];
        dDispls[i] = static_cast<int>(i ? dDispls[i - 1] + dRecvcounts[i - 1] : 0);
      }
    }
    const size_t dRecvsize = std::accumulate(dRecvcounts.begin(), dRecvcounts.end(), 0);
    distances[iDirac].resize(dRecvsize);
    values[iDirac].resize(dRecvsize);
    comm.gatherv(locDistances[iDirac], distances[iDirac], dRecvcounts, dDispls, 0);
    comm.gatherv(locValues[iDirac], values[iDirac], dRecvcounts, dDispls, 0);
  }

  // Sort by increasing separation distance
  if (comm.rank() == 0) {
    for (size_t iDirac = 0; iDirac < nDiracs; iDirac++) {
      // Get sorting order
      std::vector<size_t> indexes(distances[iDirac].size());
      std::iota(indexes.begin(), indexes.end(), 0);
      std::sort(indexes.begin(), indexes.end(),
                [&](size_t i1, size_t i2){return distances[iDirac][i1] <= distances[iDirac][i2];});

      // Apply sorting order to both distances and values
      applyIndexingDouble(indexes, distances[iDirac]);
      applyIndexingDouble(indexes, values[iDirac]);
    }
  }

  // Remove duplicates, relying on distances being ordered.
  if (removeDuplicates) {
    constexpr double relativeTolLength = 1e-5;
    const double tolLength = relativeTolLength * maxLength;
    for (size_t iDirac = 0; iDirac < nDiracs; iDirac++) {
      for (size_t node = 1; node < distances[iDirac].size(); node++) {
        if (distances[iDirac][node] - distances[iDirac][node-1] < tolLength) {
          // Values at node and node-1 are duplicates.
          // Erase value at node, and decrement node
          // so that next loop uses the new value at the same node.
          distances[iDirac].erase(distances[iDirac].begin() + node);
          values[iDirac].erase(values[iDirac].begin() + node);
          node--;
        }
      }
    }
  }

  oops::Log::trace() << "sortBySeparationDistance about to exit..." << std::endl;
  return std::tuple(distances, values, levs, fieldIndex);
}

// -----------------------------------------------------------------------------
}  // namespace util
