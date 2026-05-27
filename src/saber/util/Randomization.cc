/*
 * (C) Copyright 2026- Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/util/Randomization.h"

#include "oops/util/Logger.h"
#include "oops/util/RandomField.h"

using atlas::array::make_view;

namespace util {

// -----------------------------------------------------------------------------

void randomCtlVec(const eckit::mpi::Comm & comm,
                  const std::vector<int> & glbIndex,
                  std::vector<double> & randVecLoc) {
  oops::Log::trace() << "util::randomCtlVec starting" << std::endl;

  // Check sizes consistency
  ASSERT(glbIndex.size() == randVecLoc.size());

  // Local size
  const int nLoc = glbIndex.size();

  // Counts
  std::vector<int> counts(comm.size());
  comm.allGather(nLoc, counts.begin(), counts.end());

  // Global size
  int nGlb = 0;
  for (const auto counts : counts) {
    nGlb += counts;
  }

  // Displacements
  std::vector<int> displs(comm.size());
  displs[0] = 0;
  for (size_t jt = 0; jt < comm.size()-1; ++jt) {
    displs[jt+1] = displs[jt]+counts[jt];
  }

  // Mapping
  std::vector<int> mapping;
  if (comm.rank() == 0) {
    mapping.resize(nGlb);
  }
  comm.gatherv(glbIndex.cbegin(), glbIndex.cend(), mapping.begin(), mapping.end(),
    counts, displs, 0);

  // Generate global random vector
  std::vector<double> randVecGlb;
  if (comm.rank() == 0) {
    randVecGlb.resize(nGlb);
    util::NormalDistributionField dist(nGlb, 0.0, 1.0);
    for (int jGlb = 0; jGlb < nGlb; ++jGlb) {
      randVecGlb[jGlb] = dist[mapping[jGlb]];
    }
  }

  // Scatter random vector
  comm.scatterv(randVecGlb.cbegin(), randVecGlb.cend(), counts, displs,
    randVecLoc.begin(), randVecLoc.end(), 0);

  oops::Log::trace() << "util::randomCtlVec done" << std::endl;
}

// -----------------------------------------------------------------------------

void randomCtlVec(const eckit::mpi::Comm & comm,
                  const std::vector<int> & glbIndex,
                  atlas::Field & field) {
  oops::Log::trace() << "util::randomCtlVec starting" << std::endl;

  // Check sizes consistency
  ASSERT(glbIndex.size() == field.size());

  // Local size
  const int nLoc = glbIndex.size();

  // Create random vector
  std::vector<double> randVecLoc(nLoc);
  randomCtlVec(comm, glbIndex, randVecLoc);

  // Fill local Field
  auto view = make_view<double, 1>(field);
  for (int jLoc = 0; jLoc < nLoc; ++jLoc) {
    view(jLoc) = randVecLoc[jLoc];
  }

  oops::Log::trace() << "util::randomCtlVec done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace util

