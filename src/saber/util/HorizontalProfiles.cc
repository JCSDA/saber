/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <algorithm>

#include "atlas/array.h"
#include "atlas/util/Earth.h"

#include "eckit/geometry/Point2.h"

#include "oops/mpi/mpi.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "saber/util/HorizontalProfiles.h"

namespace util {
// -----------------------------------------------------------------------------

std::tuple<std::vector<std::vector<double>>,
           std::vector<std::vector<double>>,
           std::vector<double>,
           std::vector<double>,
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

  // Loop over fields to save (distance, value, global_index) tuples
  const auto & ghost = fspace.ghost();
  const auto ghview  = atlas::array::make_view<int, 1>(ghost);
  const auto gidxview = atlas::array::make_view<atlas::gidx_t, 1>(fspace.global_index());
  const auto lonlatView = atlas::array::make_view<double, 2>(fspace.lonlat());
  std::vector<std::vector<double>> locDistances(lons.size());
  std::vector<std::vector<double>> locValues(lons.size());
  std::vector<std::vector<atlas::gidx_t>> locGlobalIndices(lons.size());
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
              locGlobalIndices[iDirac].push_back(gidxview(jnode));
            }
          }
        }
      }
    }
  }

  // Gather distances, values, and global_indices onto root MPI task for sorting
  std::vector<std::vector<double>> distances(nDiracs);
  std::vector<std::vector<double>> values(nDiracs);
  std::vector<std::vector<atlas::gidx_t>> globalIndices(nDiracs);
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
    globalIndices[iDirac].resize(dRecvsize);
    comm.gatherv(locDistances[iDirac], distances[iDirac], dRecvcounts, dDispls, 0);
    comm.gatherv(locValues[iDirac], values[iDirac], dRecvcounts, dDispls, 0);
    comm.gatherv(locGlobalIndices[iDirac], globalIndices[iDirac], dRecvcounts, dDispls, 0);
  }

  // Sort by increasing separation distance
  if (comm.rank() == 0) {
    for (size_t iDirac = 0; iDirac < nDiracs; iDirac++) {
      // Get sorting order
      std::vector<size_t> indexes(distances[iDirac].size());
      std::iota(indexes.begin(), indexes.end(), 0);

      // Apply a fuzzy sort:
      // - if distances are more different than tolerance, sort by distance
      // - if distances are essentially equal, sort by atlas grid index
      // This makes the sort more robust against floating-point error
      std::sort(indexes.begin(), indexes.end(),
                [&](size_t i1, size_t i2) {
                  const double d1 = distances[iDirac][i1];
                  const double d2 = distances[iDirac][i2];
                  if (abs(d1 - d2) > 1e-13 * (d1 + d2)) {
                    // sufficiently different
                    return d1 <= d2;
                  } else {
                    // essentially equidistant points, fall back to grid structure
                    const atlas::gidx_t gidx1 = globalIndices[iDirac][i1];
                    const atlas::gidx_t gidx2 = globalIndices[iDirac][i2];
                    return gidx1 <= gidx2;
                  }
                });

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
  return std::tuple(distances, values, lons, lats, levs, fieldIndex);
}

// -----------------------------------------------------------------------------

void  write_1d_covariances(const eckit::mpi::Comm & comm,
                           const std::vector<std::vector<double>> & distances,
                           const std::vector<std::vector<double>> & covariances,
                           const std::vector<double> & lons,
                           const std::vector<double> & lats,
                           const std::vector<size_t> & levs,
                           const std::vector<size_t> & fieldIndexes,
                           const std::vector<std::string> & names,
                           const std::string & filePath) {
  oops::Log::trace() << "write_1d_covariances starting" << std::endl;

  // This function assumes everything has already been gathered on root PE.

  // There is no guarantee that the vectors in `distances` or `covariances` have
  // common dimensions, if they have been computed with respect to different
  // horizontal position. We need to define an independent dimension for each
  // horizontal position (lon-lat point).

  // We write the covariance profiles as 2D Fields of dimensions (horizontal bins,
  // vertical levels). In simple cases, the number of vertical levels is 1.
  // Each 2D Field is associated to a 1D variable of horizontal distances, and to
  // 1D variable of vertical level indexes.

  if (comm.rank() == 0) {
    // 1. Define variables and dimensions
    size_t nProfiles = distances.size();

    std::vector<std::string> dimNames;
    std::vector<atlas::idx_t> dimSizes;
    std::vector<std::vector<std::string>> dimNamesForEveryVar;
    std::vector<std::vector<atlas::idx_t>> dimSizesForEveryVar;
    eckit::LocalConfiguration netcdfMetaData;
    oops::Variables vars;

    std::vector<double> lonsOf2DField;
    std::vector<double> latsOf2DField;
    std::vector<std::string> varsOf2DField;

    // This vectors gives a mapping from indexes of 1D horizontal covariance profiles
    // to indexes of 2D covariance fields.
    std::vector<size_t> profilesToFieldIndex(nProfiles);

    // Absolute tolerance for longitude / latitude comparisons (in degrees)
    constexpr double lonTolerance(1e-8);
    for (size_t iProfile = 0; iProfile < nProfiles; iProfile++) {
      const std::string varName = names[fieldIndexes[iProfile]];

      // Have these horizontal position and variable already been seen or should
      // we create a new 2D field?
      bool newField(true);
      size_t iField = 0;
      while (iField < lonsOf2DField.size()) {
        if (std::abs(lonsOf2DField[iField]- lons[iProfile]) < lonTolerance) {
          if (std::abs(latsOf2DField[iField] - lats[iProfile]) < lonTolerance) {
            if (varsOf2DField[iField] == varName) {
              newField = false;
              break;
            }
          }
        }
        iField++;
      }
      profilesToFieldIndex[iProfile] = iField;

      if (newField) {
        lonsOf2DField.push_back(lons[iProfile]);
        latsOf2DField.push_back(lats[iProfile]);
        varsOf2DField.push_back(names[fieldIndexes[iProfile]]);

        // Create new horizontal dimension "ndist".
        // Dimension index in dimNames et al. will be 2*iField
        dimNames.push_back("ndist" + std::to_string(iField));
        dimSizes.push_back(distances[iProfile].size());

        // Create new vertical dimensions "levels"
        // Dimension index in dimNames et al. will be 2*iField+1
        dimNames.push_back("nz" + std::to_string(iField));
        dimSizes.push_back(1);  // May be increased later if other levels are added

        // Create new distance variable of size [ndist]
        // Variable index in vars et al. will be 3*iField
        const int varIndex = std::count(varsOf2DField.cbegin(), varsOf2DField.cend(), varName);
        const auto fieldVarSuffix = varName + "_" + std::to_string(varIndex);
        vars.push_back(fieldVarSuffix + "_distances");
        dimNamesForEveryVar.push_back({dimNames[2*iField]});
        dimSizesForEveryVar.push_back({dimSizes[2*iField]});
        util::setAttribute<std::string>(
          netcdfMetaData, vars[3*iField].name(), "binning_type", "string",
          "horizontal separation distance");
        util::setAttribute<double>(
          netcdfMetaData, vars[3*iField].name(), "longitude_of_reference_point",
          "real64", lons[iProfile]);
        util::setAttribute<double>(
          netcdfMetaData, vars[3*iField].name(), "latitude_of_reference_point",
          "real64", lats[iProfile]);

        // Create new levels variables of size [levels]
        // Variable index in vars et al. will be 3*iField+1
        vars.push_back(oops::Variable(
          fieldVarSuffix + "_levels",
          oops::VariableMetaData(oops::defaultVerticalStagger, oops::ModelDataType::Int32)));
        dimNamesForEveryVar.push_back({dimNames[2*iField+1]});
        dimSizesForEveryVar.push_back({dimSizes[2*iField+1]});

        // Create new covariance variable of size [ndist, levels]
        // Variable index in vars et al. will be 3*iField+2
        vars.push_back(fieldVarSuffix + "_covariances");
        dimNamesForEveryVar.push_back({dimNames[2*iField], dimNames[2*iField+1]});
        dimSizesForEveryVar.push_back({dimSizes[2*iField], dimSizes[2*iField+1]});
        util::setAttribute<std::string>(
          netcdfMetaData, vars[3*iField+2].name(), "statistics_type", "string",
          "1D horizontal covariance");
        util::setAttribute<double>(
          netcdfMetaData, vars[3*iField+2].name(), "longitude_of_reference_point",
          "real64", lons[iProfile]);
        util::setAttribute<double>(
          netcdfMetaData, vars[3*iField+2].name(), "latitude_of_reference_point",
          "real64", lats[iProfile]);
      } else {
        // Augment level dimensions by 1
        dimSizes[2*iField+1] += 1;
        dimSizesForEveryVar[3*iField+1][0] += 1;
        dimSizesForEveryVar[3*iField+2][1] += 1;
      }
    }

    // 2. Write Header
    std::vector<int> netcdfGeneralIDs;
    std::vector<int> netcdfDimIDs;
    std::vector<int> netcdfVarIDs;
    std::vector<std::vector<int>> netcdfDimVarIDs;
    util::atlasArrayWriteHeader(filePath,
                                dimNames,
                                dimSizes,
                                vars,
                                dimNamesForEveryVar,
                                netcdfMetaData,
                                netcdfGeneralIDs,
                                netcdfDimIDs,
                                netcdfVarIDs,
                                netcdfDimVarIDs);


    // 3. Check distance vectors associated to a given field are all equal
    const size_t nFields = dimSizes.size() / 2;
    constexpr double distAbsTol = 1e-2;  // absolute tolerance in meters for distance comparisons
    for (size_t iProfile = 0; iProfile < nProfiles; iProfile++) {
      const size_t iField = profilesToFieldIndex[iProfile];
      const auto firstProfileIndexIter = std::find(
                    profilesToFieldIndex.cbegin(),
                    profilesToFieldIndex.cbegin() + iProfile + 1,
                    iField);
      const size_t firstProfileIndex = firstProfileIndexIter - profilesToFieldIndex.cbegin();
      const bool isFirstProfile = (firstProfileIndex == iProfile);

      // Check distances are equal as we will only save one vector of distances per 2D field
      if (!isFirstProfile) {
        ASSERT(distances[iProfile].size() == distances[firstProfileIndex].size());
        for (size_t jnode = 0; jnode < distances[iProfile].size(); jnode++) {
          ASSERT(std::abs(distances[iProfile][jnode] - distances[firstProfileIndex][jnode])
                 < distAbsTol);
        }
      }
    }

    // 4. Write Data
    for (size_t iField = 0; iField< nFields; iField++) {
      const size_t ndist(dimSizes[2*iField]);
      const size_t nlevs(dimSizes[2*iField+1]);

      // Write distances
      const auto iProfileIter = std::find(profilesToFieldIndex.cbegin(),
                                          profilesToFieldIndex.cend(),
                                          iField);
      const size_t iFirstProfile = iProfileIter - profilesToFieldIndex.cbegin();
      auto distArray = atlas::array::Array::create<double>(ndist);
      auto distView = atlas::array::make_view<double, 1>(*distArray);
      for (size_t jnode = 0; jnode < ndist; jnode++) {
        distView(jnode) = distances[iFirstProfile][jnode];
      }
      auto cview = atlas::array::make_view<const double, 1>(*distArray);
      util::atlasArrayWriteData(netcdfGeneralIDs, netcdfVarIDs[3*iField], cview);

      // Create and fill levels and covariance fields
      auto levField = atlas::Field("",
                                   atlas::array::make_datatype<int>(),
                                   atlas::array::make_shape(nlevs));
      auto levView = atlas::array::make_view<int, 1>(levField);
      levView.assign(-1);

      auto covField = atlas::Field("",
                                   atlas::array::make_datatype<double>(),
                                   atlas::array::make_shape(ndist, nlevs));
      auto covView = atlas::array::make_view<double, 2>(covField);
      covView.assign(0.0);

      size_t jlev(0);
      for (size_t iProfile=0; iProfile < nProfiles; iProfile++) {
        if (profilesToFieldIndex[iProfile] == iField) {
          const int lev = levs[iProfile];
          levView(jlev) = lev;
          for (size_t jdist = 0; jdist < ndist; jdist++) {
            covView(jdist, jlev) = covariances[iProfile][jdist];
          }
          jlev++;
        }
      }

      // Write levels field
      auto constLevView = atlas::array::make_view<const int, 1>(levField);
      util::atlasArrayWriteData(netcdfGeneralIDs, netcdfVarIDs[3*iField+1], constLevView);

      // Write covariance field
      auto constCovView = atlas::array::make_view<const double, 2>(covField);
      util::atlasArrayWriteData(netcdfGeneralIDs, netcdfVarIDs[3*iField+2], constCovView);
    }

    oops::Log::info() << "Info     : covariance profiles written to file "
                      << filePath << std::endl;

    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) {
      throw eckit::Exception("NetCDF closing error", Here());
    }
  }

  oops::Log::trace() << "write_1d_covariances done" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace util
