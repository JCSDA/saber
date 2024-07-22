/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/util/BinnedFieldSetHelpers.h"

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

namespace util {

using atlas::array::make_view;
using atlas::idx_t;

// Use this version of the gatherSumFieldSet when the local bin index
// is the same as the global bin index and when the local bin index
// is the same on each PE.
atlas::FieldSet gatherSumFieldSet(const eckit::mpi::Comm & comm,
                                  const int root,
                                  const atlas::FieldSet & fset)  {
  const idx_t totalBins = static_cast<idx_t>(fset[0].shape()[0]);
  atlas::Field glBinIdx =
    atlas::Field(std::string{"binning index"}, atlas::array::make_datatype<std::int32_t>(),
                 atlas::array::make_shape(totalBins));
  auto binIndxView = make_view<std::int32_t, 1>(glBinIdx);
  for (atlas::idx_t i = 0; i < totalBins; ++i) {
    binIndxView(i) = i;
  }
  return gatherSumFieldSet(comm, root, totalBins, glBinIdx, fset);
}

// -----------------------------------------------------------------------------

atlas::FieldSet gatherSumFieldSet(const eckit::mpi::Comm & comm,
                                  const int root,
                                  const std::size_t totalbins,
                                  const atlas::Field & glBinIdx,
                                  const atlas::FieldSet & fset) {
  atlas::FieldSet gatheredFset;
  for (const atlas::Field & fld : fset) {
    if (fld.rank() == 2) {
      atlas::Field fld2 = gatherSumBinsLevels(
        comm, root, totalbins, fld.shape()[1], glBinIdx, fld);
      gatheredFset.add(fld2);
    } else if (fld.rank() == 3) {
      atlas::Field fld2 = gatherSumBinsLevelsLevels(
        comm, root, totalbins, fld.shape()[1], fld.shape()[2], glBinIdx, fld);
      gatheredFset.add(fld2);
    } else {
      oops::Log::error() << " gathering and summing a field "
                         << fld.name()
                         << " with rank "
                         << fld.rank()
                         << " is not supported." << std::endl;
    }
  }

  return gatheredFset;
}

// -----------------------------------------------------------------------------

atlas::Field gatherSumBinsLevels(
    const eckit::mpi::Comm & comm, const int root,
    const std::size_t & nbins, const std::size_t & levels,
    const atlas::Field & glBinIdx, const atlas::Field & fld) {
  const std::size_t commSize = comm.size();
  const int commRank = comm.rank();
  const int nbinsUse = (commRank == root ? nbins : 1);
  const int levelsUse = (commRank == root ? levels : 1);

  atlas::Field fldOut =
    atlas::Field(fld.name(), atlas::array::make_datatype<double>(),
                 atlas::array::make_shape(nbinsUse, levelsUse));
  //--------------------------------------------------------------------------
  // need the glBinIdx concatenated in a vector on root PE.
  //--------------------------------------------------------------------------
  auto glBinIdxView = make_view<std::int32_t, 1>(glBinIdx);
  // getting the stride and offset
  std::vector<int> glIdxStride(commSize, 0);
  std::vector<int> glIdxOffset(commSize, 0);
  std::vector<int> glIdxVec(1, glBinIdxView.size());
  int glIdxSize(0);

  comm.gather(glIdxVec, glIdxStride, root);
  if (commRank == root) {
    for (std::size_t i = 0; i < commSize; ++i) {
      glIdxOffset[i] = (i ? glIdxOffset[i - 1] + glIdxStride[i - 1] : 0);
      glIdxSize += glIdxStride[i];
    }
  }

  // collect the full size of the vector across all ranks
  std::vector<std::int32_t> glIdxRecv;
  (commRank == root ? glIdxRecv.resize(glIdxSize) : glIdxRecv.resize(0));
  std::vector<std::int32_t> glIdxSend(glBinIdxView.size(), 0);

  std::size_t n = 0;
  for (idx_t t = 0; t < glBinIdxView.shape()[0]; ++t) {
    glIdxSend[n] = glBinIdxView(t);
    ++n;
  }

  comm.gatherv(glIdxSend, glIdxRecv, glIdxStride, glIdxOffset, root);

  //--------------------------------------------------------------------------
  // need the fld gathered onto root PE.
  //--------------------------------------------------------------------------
  // collect the full size of the fld across all ranks
  auto fldView = make_view<double, 2>(fld);
  // getting the stride and offset
  std::vector<int> fldStride(commSize, 0);
  std::vector<int> fldOffset(commSize, 0);
  std::vector<int> fldVec(1, fldView.size());
  int fldSize(0);

  comm.gather(fldVec, fldStride, root);
  if (commRank == root) {
    for (std::size_t i = 0; i < commSize; ++i) {
      fldOffset[i] = (i ? fldOffset[i - 1] + fldStride[i - 1] : 0);
      fldSize += fldStride[i];
    }
  }

  std::vector<double> recv;
  std::vector<double> send(fld.size(), 0.0);
  commRank == root ? recv.resize(fldSize, 0.0) : recv.resize(0);

  n = 0;
  for (idx_t t0 = 0; t0 < fldView.shape()[0]; ++t0) {
    for (idx_t t1 = 0; t1 < fldView.shape()[1]; ++t1) {
      send[n] = fldView(t0, t1);
      ++n;
    }
  }

  comm.gatherv(send, recv, fldStride, fldOffset, root);

  auto fldOutView = make_view<double, 2>(fldOut);
  fldOutView.assign(0.0);

  if (commRank == root) {
    n = 0;
    for (std::size_t pe = 0; pe < commSize; ++pe) {
      idx_t nbinpe = static_cast<idx_t>(glIdxStride[pe]);
      for (idx_t t0 = 0; t0 < nbinpe; ++t0) {
        idx_t t2 = glIdxRecv[glIdxOffset[pe] + t0];
        for (idx_t t1 = 0; t1 < static_cast<idx_t>(levels); ++t1) {
          fldOutView(t2, t1) += recv[n];
          ++n;
        }
      }
    }
  }

  return fldOut;
}


atlas::Field gatherSumBinsLevelsLevels(
  const eckit::mpi::Comm & comm, const int root,
  const std::size_t & nbins, const std::size_t & levels1,
  const std::size_t & levels2, const atlas::Field & glBinIdx,
  const atlas::Field & fld) {
  const std::size_t commSize = comm.size();
  const int commRank = comm.rank();
  const int nbinsUse = (commRank == root ? nbins : 1);
  const int levels1Use = (commRank == root ? levels1 : 1);
  const int levels2Use = (commRank == root ? levels2 : 1);

  atlas::Field fldOut =
    atlas::Field(fld.name(), atlas::array::make_datatype<double>(),
                 atlas::array::make_shape(nbinsUse, levels1Use, levels2Use));
  //--------------------------------------------------------------------------
  // need the glBinIdx concatenated in a vector on root PE.
  //--------------------------------------------------------------------------
  auto glBinIdxView = make_view<std::int32_t, 1>(glBinIdx);
  // getting the stride and offset
  std::vector<int> glIdxStride(commSize, 0);
  std::vector<int> glIdxOffset(commSize, 0);
  std::vector<int> glIdxVec(1, glBinIdxView.size());
  int glIdxSize(0);

  comm.gather(glIdxVec, glIdxStride, root);
  if (commRank == root) {
    for (std::size_t i = 0; i < commSize; ++i) {
      glIdxOffset[i] = (i ? glIdxOffset[i - 1] + glIdxStride[i - 1] : 0);
      glIdxSize += glIdxStride[i];
    }
  }

  // collect the full size of the vector across all ranks
  std::vector<std::int32_t> glIdxRecv;
  commRank == root ? glIdxRecv.resize(glIdxSize) : glIdxRecv.resize(0);
  std::vector<std::int32_t> glIdxSend(glBinIdxView.size(), 0);

  size_t n = 0;
  for (idx_t t = 0; t < glBinIdxView.shape()[0]; ++t) {
    glIdxSend[n] = glBinIdxView(t);
    ++n;
  }

  comm.gatherv(glIdxSend, glIdxRecv, glIdxStride, glIdxOffset, root);

  //--------------------------------------------------------------------------
  // need the fld gathered onto root PE.
  //--------------------------------------------------------------------------
  // collect the full size of the fld across all ranks
  auto fldView = make_view<double, 3>(fld);
  // getting the stride and offset
  std::vector<int> fldStride(commSize, 0);
  std::vector<int> fldOffset(commSize, 0);
  std::vector<int> fldVec(1, fldView.size());
  int fldSize(0);

  comm.gather(fldVec, fldStride, root);
  if (commRank == root) {
    for (std::size_t i = 0; i < commSize; ++i) {
      fldOffset[i] = (i ? fldOffset[i - 1] + fldStride[i - 1] : 0);
      fldSize += fldStride[i];
    }
  }

  std::vector<double> recv;
  std::vector<double> send(fld.size(), 0.0);
  commRank == root ? recv.resize(fldSize, 0.0) : recv.resize(0);

  n = 0;
  for (idx_t t0 = 0; t0 < fldView.shape()[0]; ++t0) {
    for (idx_t t1 = 0; t1 < fldView.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < fldView.shape()[2]; ++t2) {
        send[n] = fldView(t0, t1, t2);
        ++n;
      }
    }
  }

  comm.gatherv(send, recv, fldStride, fldOffset, root);

  auto fldOutView = make_view<double, 3>(fldOut);
  fldOutView.assign(0.0);

  if (commRank == root) {
    n = 0;
    for (std::size_t pe = 0; pe < comm.size(); ++pe) {
      idx_t nbinpe = static_cast<idx_t>(glIdxStride[pe]);
      for (idx_t t0 = 0; t0 < nbinpe; ++t0) {
        idx_t t3 = glIdxRecv[glIdxOffset[pe] + t0];
        for (idx_t t1 = 0; t1 < static_cast<idx_t>(levels1); ++t1) {
          for (idx_t t2 = 0; t2 < static_cast<idx_t>(levels2); ++t2) {
            fldOutView(t3, t1, t2) += recv[n];
            ++n;
          }
        }
      }
    }
  }

  return fldOut;
}


}  // namespace util
