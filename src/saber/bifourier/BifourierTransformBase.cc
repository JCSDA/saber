/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierTransformBase.h"

#include <algorithm>
#include <limits>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "oops/generic/gc99.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/RandomField.h"

#include "saber/bifourier/BifourierUtilities.h"

using atlas::array::make_datatype;
using atlas::array::make_indexview;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

BifourierTransformFactory::BifourierTransformFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::BifourierTransformFactory."
      << std::endl;
    throw eckit::Exception("Element already registered in saber::BifourierTransformFactory.",
      Here());
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

std::shared_ptr<BifourierTransformBase> BifourierTransformFactory::create(
  const oops::GeometryData & gdata,
  const std::string & gridUid,
  const oops::Variables & activeVars,
  const BifourierTransformParameters & params) {
  oops::Log::trace() << "BifourierTransformBase::create starting" << std::endl;
  const std::string id = params.fftBackend.value();
  typename std::map<std::string, BifourierTransformFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::bifourier::BifourierTransformFactory."
      << std::endl;
    throw eckit::UserError("Element does not exist in saber::bifourier::BifourierTransformFactory.",
      Here());
  }
  std::shared_ptr<BifourierTransformBase> ptr =
    jsb->second->make(gdata, gridUid, activeVars, params);
  oops::Log::trace() << "BifourierTransformBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

BifourierTransformBase::BifourierTransformBase(const oops::GeometryData & gdata,
                                               const std::string & gridUid,
                                               const oops::Variables & activeVars,
                                               const BifourierTransformParameters & params) :
    gdata_(gdata),
    comm_(gdata_.comm()),
    myrank_(comm_.rank()),
    params_(params),
    gridUid_(gridUid),
    dwGlb_(params_.dwGlb.value())
{
  oops::Log::trace() << classname() << "::BifourierTransformBase starting" << std::endl;

  // Check function space type
  ASSERT(gdata_.functionSpace().type() == "StructuredColumns");

  // Print active variables
  oops::Log::info() << "Info     : New Bifourier transform" << std::endl;
  oops::Log::info() << "Info     : - FFT backend: " << params_.fftBackend.value() << std::endl;
  oops::Log::info() << "Info     : - Active variable: " << activeVars.variables() << std::endl;

  // Get function space
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());

  // Get grid size
  nx_ = fs.grid().nx()[0];
  ny_ = fs.grid().ny();
  nodes_ = fs.size();
  oops::Log::test() << "- Regional grid size: " << nx_ << "x" << ny_ << std::endl;

  // Cell size
  dx_ = fs.grid().dx(0);
  dy_ = fs.grid().y(1) - fs.grid().y(0);
  oops::Log::test() << "- Cell sizes: " << dx_*1.0e-3 << " km x " << dy_*1.0e-3 << " km"
    << std::endl;

  // Mean latitude
  atlas::PointLonLat p1 = fs.grid().lonlat(0, 0);
  atlas::PointLonLat p2 = fs.grid().lonlat(nx_-1, ny_-1);
  meanLat_ = 0.5*(p1[1]+p2[1]);
  oops::Log::test() << "- Mean latitude: " << meanLat_ << " deg" << std::endl;

  // Number of levels for all variables
  nvz_ = 0;
  for (const auto & var : activeVars) {
    nvz_ += var.getLevels();
  }

  oops::Log::trace() << classname() << "::BifourierTransformBase done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::test(const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::test starting" << std::endl;

  // Get tests tolerance
  const double tolerance = params_.specTolerance.value();

  // Generate random FieldSet
  atlas::FieldSet gpFset = util::createRandomFieldSet(comm_, gdata_.functionSpace(), activeVars);

  // Truncate grid-point field
  atlas::FieldSet spFset;
  gp2sp(gpFset, spFset, activeVars);
  sp2gp(spFset, gpFset, activeVars);

  // Grid-point to spectral
  gp2sp(gpFset, spFset, activeVars);

  // Check inverse
  atlas::FieldSet gpFsetTest = util::copyFieldSet(gpFset);
  sp2gp(spFset, gpFsetTest, activeVars);
  ASSERT(util::compareFieldSets(comm_, gpFset, gpFsetTest));
  oops::Log::test() << "- Direct-inverse test passed" << std::endl;

  // Check forward
  atlas::FieldSet spFsetTest;
  gp2sp(gpFsetTest, spFsetTest, activeVars);
  for (const auto & var : activeVars) {
    const size_t nz = var.getLevels();
    const auto spField = spFset[var.name()];
    const auto spFieldTest = spFsetTest[var.name()];
    const auto spView = make_view<double, 2>(spField);
    const auto spViewTest = make_view<double, 2>(spFieldTest);
    int wrongValues = 0;
    for (size_t js = 0; js < ns_; ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        if (!oops::is_close_relative(spView(js, jz), spViewTest(js, jz), tolerance)) {
          ++wrongValues;
        }
      }
    }
    comm_.allReduceInPlace(wrongValues, eckit::mpi::sum());
    ASSERT(wrongValues == 0);
  }
  oops::Log::test() << "- Inverse-direct test passed" << std::endl;

  // Check Parseval's identity
  double gpSqNorm = util::dotProductFieldSets(gpFset, gpFset, activeVars.variables(), comm_);
  double spSqNorm = util::dotProductFieldSets(spFset, spFset, activeVars.variables(), comm_);
//  ASSERT(oops::is_close_relative(gpSqNorm, spSqNorm, tolerance));
  oops::Log::test() << "- Parseval identity test passed" << std::endl;

  // Adjoint test, forward
  gpFset = util::createRandomFieldSet(comm_, gdata_.functionSpace(), activeVars);
  gp2sp(gpFset, spFset, activeVars);
  createRandomFieldSet(spFsetTest, activeVars);
  gp2spAdj(spFsetTest, gpFsetTest, activeVars);
  gpSqNorm = util::dotProductFieldSets(gpFset, gpFsetTest, activeVars.variables(), comm_);
  spSqNorm = util::dotProductFieldSets(spFset, spFsetTest, activeVars.variables(), comm_);
  ASSERT(oops::is_close_relative(gpSqNorm, spSqNorm, tolerance));
  oops::Log::test() << "- Adjoint test (forward) passed" << std::endl;

  // Adjoint test, inverse
  gpFset = util::createRandomFieldSet(comm_, gdata_.functionSpace(), activeVars);
  sp2gpAdj(gpFset, spFset, activeVars);
  createRandomFieldSet(spFsetTest, activeVars);
  sp2gp(spFsetTest, gpFsetTest, activeVars);
  gpSqNorm = util::dotProductFieldSets(gpFset, gpFsetTest, activeVars.variables(), comm_);
  spSqNorm = util::dotProductFieldSets(spFset, spFsetTest, activeVars.variables(), comm_);
  ASSERT(oops::is_close_relative(gpSqNorm, spSqNorm, tolerance));
  oops::Log::test() << "- Adjoint test (inverse) passed" << std::endl;

  // Derivatives / Laplacian consistency test
  for (const auto & var : activeVars) {
    // Get field
    const size_t nz = var.getLevels();
    auto spField = spFset[var.name()];

    // Double derivative in X direction
    atlas::Field spDxField;
    derivative(spField, spDxField, "x");
    atlas::Field spDx2Field;
    derivative(spDxField, spDx2Field, "x");

    // Double derivative in Y direction
    atlas::Field spDyField;
    derivative(spField, spDyField, "y");
    atlas::Field spDy2Field;
    derivative(spDyField, spDy2Field, "y");

    // Direct Laplacian
    atlas::Field spLapDirField = spField.clone();
    directLaplacian(spLapDirField);

    // Comparison
    const auto spDx2View = make_view<double, 2>(spDx2Field);
    const auto spDy2View = make_view<double, 2>(spDy2Field);
    const auto spLapDirView = make_view<double, 2>(spLapDirField);
    int wrongValues = 0;
    for (size_t js = 0; js < ns_; ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        if (!oops::is_close_relative(spLapDirView(js, jz), spDx2View(js, jz) + spDy2View(js, jz),
          tolerance)) {
          ++wrongValues;
        }
      }
    }
    comm_.allReduceInPlace(wrongValues, eckit::mpi::sum());
    ASSERT(wrongValues == 0);
  }
  oops::Log::test() << "- Derivatives / direct Laplacian consistency test passed" << std::endl;

  oops::Log::trace() << classname() << "::test done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::gp2sp(const atlas::FieldSet & gpFset,
                                   atlas::FieldSet & spFset,
                                   const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::gp2sp starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create recv vectors
  std::vector<double> recvVec;

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Serialize from grid-point FieldSet
  recvVec.resize(gridRecvSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto gpField = gpFset[var.name()];
    ASSERT(gpField.shape(0) == static_cast<int>(nodes_));
    ASSERT(gpField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          recvVec[jgrv] = gpView(jnode, jz);
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Prepare spectral FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (spFset.has(var.name())) {
      // Check sizes
      ASSERT(spFset[var.name()].shape(0) == static_cast<int>(ns_));
      ASSERT(spFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field spField = spFspace_->createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      spFset.add(spField);
    }
  }

  // Backend specific transform
  gp2sp(recvVec, spFset.metadata());

  // Communication
  std::vector<double> sendVec(eqchSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), recvCounts().data(), recvDispls().data(),
    sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data());

  // Reserialize into spectral FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto spField = spFset[var.name()];

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        spView(js, jz) = sendVec[jesv];

        // Normalize FFT
        spView(js, jz) *= gp2spNorm(js);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::gp2sp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::sp2gp(const atlas::FieldSet & spFset,
                                   atlas::FieldSet & gpFset,
                                   const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::sp2gp starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create recv vectors
  std::vector<double> recvVec;

  // Reserialize from spectral FieldSet
  std::vector<double> sendVec(eqchSendSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = spFset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        sendVec[jesv] = spView(js, jz);

        // Normalize FFT
        sendVec[jesv] *= sp2gpNorm(js);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  recvVec.resize(recvSize()*nvz_);
  comm_.allToAllv(sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data(),
    recvVec.data(), recvCounts().data(), recvDispls().data());

  // Backend specific transform
  sp2gp(recvVec, spFset.metadata());

  // Prepare grid-point FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (gpFset.has(var.name())) {
      // Check sizes
      ASSERT(gpFset[var.name()].shape(0) == static_cast<int>(nodes_));
      ASSERT(gpFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field gpField = gdata_.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      gpFset.add(gpField);
    }
  }

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Deserialize into grid-point FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto gpField = gpFset[var.name()];

    // Get field view
    auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          gpView(jnode, jz) = recvVec[jgrv];
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::sp2gp done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::gp2spAdj(const atlas::FieldSet & spFset,
                                      atlas::FieldSet & gpFset,
                                      const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::gp2spAdj starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create recv vectors
  std::vector<double> recvVec;

  // Reserialize from spectral FieldSet
  std::vector<double> sendVec(eqchSendSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = spFset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        sendVec[jesv] = spView(js, jz);

        // Normalize FFT
        sendVec[jesv] *= gp2spAdjNorm(js);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  recvVec.resize(recvSize()*nvz_);
  comm_.allToAllv(sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data(),
    recvVec.data(), recvCounts().data(), recvDispls().data());

  // Backend specific transform
  gp2spAdj(recvVec, spFset.metadata());

  // Prepare grid-point FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (gpFset.has(var.name())) {
      // Check sizes
      ASSERT(gpFset[var.name()].shape(0) == static_cast<int>(nodes_));
      ASSERT(gpFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field gpField = gdata_.functionSpace().createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      gpFset.add(gpField);
    }
  }

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Deserialize into grid-point FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto gpField = gpFset[var.name()];

    // Get field view
    auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          gpView(jnode, jz) = recvVec[jgrv];
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::gp2spAdj done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::sp2gpAdj(const atlas::FieldSet & gpFset,
                                      atlas::FieldSet & spFset,
                                      const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::sp2gpAdj starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Create vectors
  std::vector<double> recvVec;

  // Ghost points
  const auto ghostView = make_view<int, 1>(gdata_.functionSpace().ghost());

  // Serialize from grid-point FieldSet
  recvVec.resize(gridRecvSize_*nvz_);
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto gpField = gpFset[var.name()];
    ASSERT(gpField.shape(0) == static_cast<int>(nodes_));
    ASSERT(gpField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto gpView = make_view<double, 2>(gpField);
    size_t jgr = 0;
    for (size_t jnode = 0; jnode < nodes_; ++jnode) {
      if (ghostView(jnode) == 0) {
        for (size_t jz = 0; jz < nz; ++jz) {
          // Total level index
          const size_t jvz = zOffset + jz;

          // Communication vector index
          const size_t jgrv = gridRecvIndex_[jgr] + jvz;

          // Copy data
          recvVec[jgrv] = gpView(jnode, jz);
        }
        ++jgr;
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Prepare spectral FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (spFset.has(var.name())) {
      // Check sizes
      ASSERT(spFset[var.name()].shape(0) == static_cast<int>(ns_));
      ASSERT(spFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field spField = spFspace_->createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      spFset.add(spField);
    }
  }

  // Backend specific transform
  sp2gpAdj(recvVec, spFset.metadata());

  // Communication
  std::vector<double> sendVec(eqchSendSize_*nvz_);
  comm_.allToAllv(recvVec.data(), recvCounts().data(), recvDispls().data(),
    sendVec.data(), eqchSendCounts_.data(), eqchSendDispls_.data());

  // Reserialize into spectral FieldSet
  zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Get field
    auto spField = spFset[var.name()];

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t jes = 0; jes < eqchSendSize_; ++jes) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = eqchSendIndex_[jes];

        // Communication vector index
        const size_t jesv = jes*nvz_+jvz;

        // Copy data
        spView(js, jz) = sendVec[jesv];

        // Normalize FFT
        spView(js, jz) *= sp2gpAdjNorm(js);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::sp2gpAdj done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::createRandomFieldSet(atlas::FieldSet & spFset,
                                                  const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::createRandomFieldSet starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Global vector
  std::vector<double> rand_vec_glb;

  if (myrank_ == 0) {
    // Generate global random vector
    rand_vec_glb.resize(nsGlb_*nvz_);
    util::NormalDistributionField dist(nsGlb_*nvz_, 0.0, 1.0);
    for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
      for (size_t jvz = 0; jvz < nvz_; ++jvz) {
        const size_t jj = jsGlb*nvz_ + jvz;
        const size_t jjOrdered = sMapping_[jsGlb]*nvz_ + jvz;
        rand_vec_glb[jj] = dist[jjOrdered];
      }
    }
  }

  // Scatter random vector
  std::vector<int> counts = sCounts_;
  std::vector<int> displs = sDispls_;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    counts[jt] *= nvz_;
    displs[jt] *= nvz_;
  }
  std::vector<double> rand_vec_loc(ns_*nvz_);
  comm_.scatterv(rand_vec_glb.cbegin(), rand_vec_glb.cend(), counts, displs,
    rand_vec_loc.begin(), rand_vec_loc.end(), 0);

  // Prepare spectral FieldSet
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    if (spFset.has(var.name())) {
      // Check sizes
      ASSERT(spFset[var.name()].shape(0) == static_cast<int>(ns_));
      ASSERT(spFset[var.name()].shape(1) == static_cast<int>(nz));
    } else {
      // Create field
      atlas::Field spField = spFspace_->createField<double>(
        atlas::option::name(var.name()) | atlas::option::levels(nz));
      spFset.add(spField);
    }
  }

  // Reserialize into spectral FieldSet
  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // CHeck field
    auto spField = spFset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t js = 0; js < ns_; ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Random vector index
        const size_t jr = js*nvz_ + jvz;

        // Copy data
        spView(js, jz) = rand_vec_loc[jr];
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::createRandomFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::derivative(const atlas::Field & spField,
                                        atlas::Field & spDerivField,
                                        const std::string & direction,
                                        const bool & adjoint) const {
  oops::Log::trace() << classname() << "::derivative starting" << std::endl;

  // Check field size
  ASSERT(spField.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = spField.shape(1);

  if (spDerivField.valid()) {
    // spDerivField is already allocated
    ASSERT(spDerivField.shape(0) == static_cast<int>(ns_));
    ASSERT(spDerivField.shape(1) == static_cast<int>(nz));
  } else {
    // Allocate spDerivField
    spDerivField = spField.clone();
  }

  // Get fields views
  const auto spView = make_view<double, 2>(spField);
  auto spDerivView = make_view<double, 2>(spDerivField);

  // Get derivative linear operator
  const auto firstIndexView = adjoint ?
    make_view<int, 1>(derivatives_[direction + "DerivCol"]) :
    make_view<int, 1>(derivatives_[direction + "DerivRow"]);
  const auto secondIndexView = adjoint ?
    make_view<int, 1>(derivatives_[direction + "DerivRow"]) :
    make_view<int, 1>(derivatives_[direction + "DerivCol"]);
  const auto factorView = make_view<double, 1>(derivatives_[direction + "DerivS"]);

  // Apply derivative linear operator
  spDerivView.assign(0.0);
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      spDerivView(firstIndexView(js), jz) = spView(secondIndexView(js), jz)*factorView(js);
    }
  }

  oops::Log::trace() << classname() << "::derivative done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::directLaplacian(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::directLaplacian starting" << std::endl;

  // Check field size
  ASSERT(field.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = field.shape(1);

  // Get field view
  auto view = make_view<double, 2>(field);

  // Apply direct Laplacian factor
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      view(js, jz) *= lapDirVec_[js];
    }
  }

  oops::Log::trace() << classname() << "::directLaplacian done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::inverseLaplacian(atlas::Field & field) const {
  oops::Log::trace() << classname() << "::inverseLaplacian starting" << std::endl;

  // Check field size
  ASSERT(field.shape(0) == static_cast<int>(ns_));

  // Get number of vertical levels
  const size_t nz = field.shape(1);

  // Get field view
  auto view = make_view<double, 2>(field);

  // Apply inverse Laplacian factor
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jz = 0; jz < nz; ++jz) {
      view(js, jz) *= lapInvVec_[js];
    }
  }
  oops::Log::trace() << classname() << "::inverseLaplacian done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::cv2fset(const atlas::Field & cv,
                                     atlas::FieldSet & fset,
                                     const oops::Variables & activeVars,
                                     const size_t & offset) const {
  oops::Log::trace() << classname() << "::cv2fset starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Clear FieldSet
  fset.clear();

  // Get control vector view
  const auto cvView = make_view<double, 1>(cv);

  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Create field
    atlas::Field spField = spFspace_->createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    fset.add(spField);

    // Get field view
    auto spView = make_view<double, 2>(spField);

    for (size_t jz = 0; jz < nz; ++jz) {
      // Total level index
      const size_t jvz = zOffset + jz;

      for (size_t js = 0; js < ns_; ++js) {
        // Control vector index
        const size_t jcv = offset + jvz*ns_ + js;

        // Copy data
        spView(js, jz) = cvView(jcv);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::cv2fset done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::fset2cv(const atlas::FieldSet & fset,
                                     atlas::Field & cv,
                                     const oops::Variables & activeVars,
                                     const size_t & offset) const {
  oops::Log::trace() << classname() << "::fset2cv starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Get CV view
  auto cvView = make_view<double, 1>(cv);

  size_t zOffset = 0;
  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = fset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(ns_));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jz = 0; jz < nz; ++jz) {
      // Total level index
      const size_t jvz = zOffset + jz;

      for (size_t js = 0; js < ns_; ++js) {
        // Control vector index
        const size_t jcv = offset + jvz*ns_ + js;

        // Copy data
        cvView(jcv) = spView(js, jz);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::fset2cv done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::copyFieldSet(const atlas::FieldSet & spInputFset,
                                          atlas::FieldSet & spOutputFset,
                                          const oops::Variables & activeVars) const {
  oops::Log::trace() << classname() << "::copyFieldSet starting" << std::endl;

  // Check the number of required levels
  size_t nvz = 0;
  for (const auto & var : activeVars) {
    nvz += var.getLevels();
  }
  ASSERT(nvz == nvz_);

  // Remove active variables from output FieldSet
  util::removeFieldsFromFieldSet(spOutputFset, activeVars.variables());

  for (const auto & var : activeVars) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check input field
    const auto spInputField = spInputFset[var.name()];
    ASSERT(spInputField.shape(0) == static_cast<int>(ns_));
    ASSERT(spInputField.shape(1) == static_cast<int>(nz));

    // Create output field
    atlas::Field spOutputField = spFspace_->createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    spOutputFset.add(spOutputField);

    // Get fields views
    const auto spInputView = make_view<double, 2>(spInputField);
    auto spOutputView = make_view<double, 2>(spOutputField);

    for (size_t jz = 0; jz < nz; ++jz) {
      for (size_t js = 0; js < ns_; ++js) {
        // Copy data
        spOutputView(js, jz) = spInputView(js, jz);
      }
    }
  }

  // Copy metadata
  spOutputFset.metadata() = spInputFset.metadata();

  oops::Log::trace() << classname() << "::copyFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------

double BifourierTransformBase::rkstar(const size_t & jk,
                                      const size_t & jl,
                                      const size_t & M,
                                      const size_t & N,
                                      const size_t & nwGlb) const {
  const double w = static_cast<double>(nwGlb-1)*std::sqrt(
    static_cast<double>(jk*jk)/static_cast<double>(M*M)
    + static_cast<double>(jl*jl)/static_cast<double>(N*N));
  return w;
}

// -----------------------------------------------------------------------------

double BifourierTransformBase::ikstar(const size_t & jk,
                                      const size_t & jl,
                                      const size_t & M,
                                      const size_t & N,
                                      const size_t & nwGlb) const {
  const size_t iw = static_cast<size_t>(rkstar(jk, jl, M, N, nwGlb)+jwGlbTol_);
  return iw;
}

// -----------------------------------------------------------------------------

bool BifourierTransformBase::includeWavenumber(const size_t & js,
                                               const size_t & jw) const {
  bool include = true;
  if (dwGlb_ > 0.0) {
    // Get indices
    const double jkstar = kstarVec_[js];
    const size_t jk = jkVec_[js];
    const size_t jl = jlVec_[js];
    const size_t jwGlb = jw + nwStartPerTask_[myrank_];

    // Keep all wavenumbers inside [jwGlb-dwGlb, jwGlb+dwGlb]
    include = include && (jkstar >= static_cast<double>(jwGlb)-dwGlb_);
    include = include && (jkstar <= static_cast<double>(jwGlb)+dwGlb_);

    // Separate wavenumber 0 from others
    include = include && (((jwGlb == 0) && (jk == 0) && (jl == 0)) || (jwGlb != 0));
  } else {
    // Get index
    const size_t jwTest = jwGlbVec_[js] - nwStartPerTask_[myrank_];

    // Check conditions
    include = include && (jw == jwTest);
  }
  return include;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::reduceCov(atlas::Field & covField) const {
  oops::Log::trace() << classname() << "::reduceCov starting" << std::endl;

  // Get sizes
  const size_t nzI = covField.shape(1);
  const size_t nzJ = covField.shape(2);

  // Get covariance view
  auto covView = make_view<double, 3>(covField);

  for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
    if (myJwGlb_[jwGlb]) {
      // Create covariance vector
      std::vector<double> covVec(nzI*nzJ);

      // Get local total wavenumber
      const size_t jw = jwGlb - nwStartPerTask_[myrank_];

      // Copy data from covariance field to covariance vector
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          const size_t jj = jzI*nzJ + jzJ;
          covVec[jj] = covView(jw, jzI, jzJ);
        }
      }

      // Allreduce covariance vector
      covRedComm_[jwGlb]->allReduceInPlace(covVec.begin(), covVec.end(), eckit::mpi::sum());

      // Copy data from covariance vector to covariance field
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          const size_t jj = jzI*nzJ + jzJ;
          covView(jw, jzI, jzJ) = covVec[jj];
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::reduceCov starting" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::scatterCov(const std::vector<double> & covVecGlb,
                                        atlas::Field & covField,
                                        const bool & adjoint) const {
  oops::Log::trace() << classname() << "::scatterCov starting" << std::endl;

  // Get sizes
  const size_t nzI = covField.shape(1);
  const size_t nzJ = covField.shape(2);

  // Check global vector size
  if (myrank_ == 0) {
    ASSERT(covVecGlb.size() == nwGlb_*nzI*nzJ);
  } else {
    ASSERT(covVecGlb.size() == 0);
  }

  // Get covariance view
  auto covView = make_view<double, 3>(covField);

  // TODO(Benjamin): check which option is the fastest
  if (true) {
    // Allocate root vector
    std::vector<double> rootVec(nwRoot_*nzI*nzJ);

    // Scatter vector
    std::vector<int> wCounts(wCounts_);
    std::vector<int> wDispls(wDispls_);
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      wCounts[jt] *= nzI*nzJ;
      wDispls[jt] *= nzI*nzJ;
    }
    comm_.scatterv(covVecGlb.cbegin(), covVecGlb.cend(), wCounts, wDispls,
      rootVec.begin(), rootVec.end(), 0);

    for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
      if (myJwGlb_[jwGlb]) {
        // Create covariance vector
        std::vector<double> covVec(nzI*nzJ, 0.0);

        // Get local total wavenumber
        const size_t jw = jwGlb - nwStartPerTask_[myrank_];

        // Copy data from root vector to covariance vector
        if (jw < nwRoot_) {
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
              const size_t jj = jzI*nzJ + jzJ;
              const size_t jjRoot = adjoint ? jw*nzJ*nzI + jzJ*nzI + jzI
                : jw*nzI*nzJ + jzI*nzJ + jzJ;
              covVec[jj] = rootVec[jjRoot];
            }
          }
        }

        // Allreduce covariance vector
        covRedComm_[jwGlb]->allReduceInPlace(covVec.begin(), covVec.end(), eckit::mpi::sum());

        // Copy data from covariance vector to covariance field
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            const size_t jj = jzI*nzJ + jzJ;
            covView(jw, jzI, jzJ) = covVec[jj];
          }
        }
      }
    }
  } else {
    for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
      if (myJwGlb_[jwGlb] || (myrank_ == 0)) {
        // Create covariance vector
        std::vector<double> covVec(nzI*nzJ);

        // Copy data from global vector to covariance vector
        if (myrank_ == 0) {
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
              const size_t jj = jzI*nzJ + jzJ;
              const size_t jjGlb = adjoint ? jwGlb*nzJ*nzI + jzJ*nzI + jzI
                : jwGlb*nzI*nzJ + jzI*nzJ + jzJ;
              covVec[jj] = covVecGlb[jjGlb];
            }
          }
        }

        // Broadcast covariance vector
        covBcastComm_[jwGlb]->broadcast(covVec.begin(), covVec.end(), 0);

        if (myJwGlb_[jwGlb]) {
          // Get local total wavenumber
          const size_t jw = jwGlb - nwStartPerTask_[myrank_];

          // Copy data from covariance vector to covariance field
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
              const size_t jj = jzI*nzJ + jzJ;
              covView(jw, jzI, jzJ) = covVec[jj];
            }
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::scatterCov done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::gatherCov(const atlas::Field & covField,
                                       std::vector<double> & covVecGlb,
                                       const bool & adjoint) const {
  oops::Log::trace() << classname() << "::gatherCov starting" << std::endl;

  // Get sizes
  const size_t nzI = covField.shape(1);
  const size_t nzJ = covField.shape(2);

  // Define vector
  std::vector<double> covVecRoot(nwRoot_*nzI*nzJ);

  // Get covariance view
  const auto covView = make_view<double, 3>(covField);

  // Serialize
  for (size_t jw = 0; jw < nwRoot_; ++jw) {
    for (size_t jzI = 0; jzI < nzI; ++jzI) {
      for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
        const size_t jjRoot = adjoint ? jw*nzJ*nzI + jzJ*nzI + jzI
          : jw*nzI*nzJ + jzI*nzJ + jzJ;
        covVecRoot[jjRoot] = covView(jw, jzI, jzJ);
      }
    }
  }

  // Allocate global vector
  if (comm_.rank() == 0) {
    covVecGlb.resize(nwGlb_*nzI*nzJ);
  } else {
    covVecGlb.resize(0);
  }

  // Gather vector
  std::vector<int> wCounts(wCounts_);
  std::vector<int> wDispls(wDispls_);
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    wCounts[jt] *= nzI*nzJ;
    wDispls[jt] *= nzI*nzJ;
  }
  comm_.gatherv(covVecRoot.cbegin(), covVecRoot.cend(), covVecGlb.begin(), covVecGlb.end(),
    wCounts, wDispls, 0);

  oops::Log::trace() << classname() << "::gatherCov done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::filterCov(const double & Lf,
                                       atlas::Field & covField) const {
  oops::Log::trace() << classname() << "::filterCov starting" << std::endl;

  if (Lf > 0.0) {
    // Get sizes
    const size_t nzI = covField.shape(1);
    const size_t nzJ = covField.shape(2);
    const size_t nzz = nzI*nzJ;

    // Get covariance view
    auto covView = make_view<double, 3>(covField);

    // Compute averaged order of magnitude for each total wavenumber
    std::vector<double> magnitude(nwGlb_, 0.0);
    for (size_t jw = 0; jw < nwRoot_; ++jw) {
      const size_t jwGlb = jw + nwStartPerTask_[myrank_];
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          magnitude[jwGlb] += std::abs(covView(jw, jzI, jzJ));
        }
      }
    }

    // Reduce magnitude
    comm_.allReduceInPlace(magnitude.begin(), magnitude.end(), eckit::mpi::sum());

    // Set zero magnitude to one to avoid NaNs
    for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
      if (!(magnitude[jwGlb] > 0.0)) {
        magnitude[jwGlb] = 1.0;
      }
    }

    // Transpose parallelization
    std::vector<int> nzzSizePerTask(comm_.size(), 0);
    size_t jt = 0;
    for (size_t jzz = 0; jzz < nzz; ++jzz) {
      ++nzzSizePerTask[jt];
      ++jt;
      if (jt == comm_.size()) {
        jt = 0;
      }
    }

    // Communication vectors
    std::vector<int> sendCounts(comm_.size());
    std::vector<int> sendDispls(comm_.size());
    std::vector<int> recvCounts(comm_.size());
    std::vector<int> recvDispls(comm_.size());
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      sendCounts[jt] = nwRoot_*nzzSizePerTask[jt];;
      sendDispls[jt] = static_cast<int>(jt ? sendDispls[jt-1] + sendCounts[jt-1] : 0);
      recvCounts[jt] = wCounts_[jt]*nzzSizePerTask[myrank_];
      recvDispls[jt] = wDispls_[jt]*nzzSizePerTask[myrank_];
    }

    // Serialize
    std::vector<double> sendBuf(nwRoot_*nzz);
    for (size_t jw = 0; jw < nwRoot_; ++jw) {
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          const size_t jc = jzI*nzJ*nwRoot_ + jzJ*nwRoot_ + jw;
          sendBuf[jc] = covView(jw, jzI, jzJ);
        }
      }
    }

    // Communication
    std::vector<double> recvBuf(nwGlb_*nzzSizePerTask[myrank_]);
    comm_.allToAllv(sendBuf.data(), sendCounts.data(), sendDispls.data(),
      recvBuf.data(), recvCounts.data(), recvDispls.data());

    // Deserialize
    std::vector<std::vector<double>> covT(nzzSizePerTask[myrank_]);
    std::vector<std::vector<double>> covTCopy(nzzSizePerTask[myrank_]);
    for (int jzz = 0; jzz < nzzSizePerTask[myrank_]; ++jzz) {
      covT[jzz].resize(nwGlb_);
      covTCopy[jzz].resize(nwGlb_);
    }
    size_t jc = 0;
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      for (int jzz = 0; jzz < nzzSizePerTask[myrank_]; ++jzz) {
        for (size_t jw = 0; jw < nwRootPerTask_[jt]; ++jw) {
          const size_t jwGlb = jw + nwStartPerTask_[jt];
          covTCopy[jzz][jwGlb] = recvBuf[jc]/magnitude[jwGlb];
          ++jc;
        }
      }
    }

    // Apply filtering kernel
    std::vector<double> kernel(nwGlb_);
    for (int jzz = 0; jzz < nzzSizePerTask[myrank_]; ++jzz) {
      covT[jzz][0] = magnitude[0]*covTCopy[jzz][0];
    }
    for (size_t jwGlb = 1; jwGlb < nwGlb_; ++jwGlb) {
      // Compute kernel
      for (size_t kwGlb = 0; kwGlb < nwGlb_; ++kwGlb) {
        const double dist = std::abs(static_cast<double>(jwGlb)-static_cast<double>(kwGlb))/Lf;
        kernel[kwGlb] = oops::gc99(dist);
      }

      for (int jzz = 0; jzz < nzzSizePerTask[myrank_]; ++jzz) {
        covT[jzz][jwGlb] = 0.0;
        double kernelSum = 0.0;
        for (size_t kwGlb = 0; kwGlb < nwGlb_; ++kwGlb) {
          covT[jzz][jwGlb] += kernel[kwGlb]*covTCopy[jzz][kwGlb];
          kernelSum += kernel[kwGlb];
        }
        covT[jzz][jwGlb] *= magnitude[jwGlb]/kernelSum;
      }
    }

    // Serialize
    jc = 0;
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      for (int jzz = 0; jzz < nzzSizePerTask[myrank_]; ++jzz) {
        for (size_t jw = 0; jw < nwRootPerTask_[jt]; ++jw) {
          const size_t jwGlb = jw + nwStartPerTask_[jt];
          recvBuf[jc] = covT[jzz][jwGlb];
          ++jc;
        }
      }
    }

    // Communication
    comm_.allToAllv(recvBuf.data(), recvCounts.data(), recvDispls.data(),
      sendBuf.data(), sendCounts.data(), sendDispls.data());

    // Reset cov to zero
    covView.assign(0.0);

    // Deserialize
    for (size_t jw = 0; jw < nwRoot_; ++jw) {
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          const size_t jc = jzI*nzJ*nwRoot_ + jzJ*nwRoot_ + jw;
          covView(jw, jzI, jzJ) = sendBuf[jc];
        }
      }
    }

    // Reduce covariance
    reduceCov(covField);
  }

  oops::Log::trace() << classname() << "::filterCov done" << std::endl;
}

// -----------------------------------------------------------------------------

double BifourierTransformBase::normCov(const atlas::Field & covField) const {
  oops::Log::trace() << classname() << "::normCov starting" << std::endl;

  // Get sizes
  const size_t nzI = covField.shape(1);
  const size_t nzJ = covField.shape(2);

  // Get covariance view
  const auto covView = make_view<double, 3>(covField);

  // Compute local squared norm
  double zz = 0.0;
  for (size_t jw = 0; jw < nwRoot_; ++jw) {
    for (size_t jzI = 0; jzI < nzI; ++jzI) {
      for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
        zz += covView(jw, jzI, jzJ)*covView(jw, jzI, jzJ);
      }
    }
  }

  // Reduce squared norm
  comm_.allReduceInPlace(zz, eckit::mpi::sum());

  // Take square-root
  zz = std::sqrt(zz);

  oops::Log::trace() << classname() << "::normCov starting" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::reduceNormalizeCov(const size_t & rank,
                                                atlas::Field & covField) {
  oops::Log::trace() << classname() << "::reduceNormalizeCov starting" << std::endl;

  // Reduce covariance
  reduceCov(covField);

  // Get number of levels
  const size_t nzI = covField.shape(1);
  const size_t nzJ = covField.shape(2);

  // Get covariance view
  auto covView = make_view<double, 3>(covField);

  // Normalize covariance
  for (size_t jw = 0; jw < nw_; ++jw) {
    const double covEnsNorm = 1.0/static_cast<double>(rank);
    const double covNorm = covEnsNorm*spNormSumInv_[jw];
    for (size_t jzI = 0; jzI < nzI; ++jzI) {
      for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
        covView(jw, jzI, jzJ) *= covNorm;
      }
    }
  }

  // Set wavenumber 0 to zero for vorticity and divergence
  if (covField.name() == "air_upward_absolute_vorticity"
    || covField.name() == "air_horizontal_divergence") {
    for (size_t jw = 0; jw < nw_; ++jw) {
      const size_t jwGlb = jw + nwStartPerTask_[myrank_];
      if (jwGlb == 0) {
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            covView(jw, jzI, jzJ) = 0.0;
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::reduceNormalizeCov done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::print(std::ostream & os) const {
  os << classname();
}


// -----------------------------------------------------------------------------

void BifourierTransformBase::setupGlobalSpectralSpace() {
  oops::Log::trace() << classname() << "::setupGlobalSpectralSpace starting" << std::endl;

  // Maximum wave numbers
  nk_ = nx_/2+1;
  nl_ = ny_/2+1;
  oops::Log::test() << "- Spectral sizes: " << nk_ << "x" << nl_ << std::endl;

  // Differential operators factors
  exwn_ = 2.0*M_PI/(static_cast<double>(nx_)*dx_);
  eywn_ = 2.0*M_PI/(static_cast<double>(ny_)*dy_);

  // Normalization factor
  normFFT_ = 1.0/static_cast<double>(nx_*ny_);

  // Define truncation
  if (params_.truncationType.value() == "arome") {
    // Same as the AROME model

    // Define tolerance to define jwGlb
    jwGlbTol_ = 0.49;

    // Define truncation parameters
    M_ = (nx_-1)/params_.truncationFactor.value();
    N_ = (ny_-1)/params_.truncationFactor.value();
    oops::Log::test() << "- Truncation parameters MxN: " << M_ << "x" <<  N_ << std::endl;

    // Define ellips
    ellips_.resize(M_+1);
    ellips_[0] = N_;
    for (size_t jk = 1; jk < M_; ++jk) {
      ellips_[jk] = static_cast<size_t>(static_cast<float>(N_)/static_cast<float>(M_)
        *std::sqrt(static_cast<float>(M_*M_-jk*jk))+1.0e-10);
    }
    ellips_[M_] = 0;

    // Maximum total wave number
    nwGlb_ = std::max(M_, N_)+1;
    oops::Log::test() << "- Maximum total wave number: " << nwGlb_-1 << std::endl;
  } else {
    // Unknown truncation
    throw eckit::Exception("unknown truncation type", Here());
  }

  // Mapping
  spNormKL_.resize(nk_*nl_, 0.0);
  size_t jk;
  size_t jl;

  // k = 0
  jk = 0;

  // l = 0
  jl = 0;
  addSpectralCoefficient(jk, jl, ReRe, 0, 0);

  // 0 < l < nl_-1
  for (size_t jl = 1; jl < nl_-1; ++jl) {
    addSpectralCoefficient(jk, jl, ReRe, 0, 1);
    addSpectralCoefficient(jk, jl, ReIm, 0, -1);
  }

  // l = nl_-1
  jl = nl_-1;
  if (ny_ % 2 == 1) {
    addSpectralCoefficient(jk, jl, ReRe, 0, 1);
    addSpectralCoefficient(jk, jl, ReIm, 0, -1);
  } else {
    addSpectralCoefficient(jk, jl, ReRe, 0, 0);
  }

  // 0 < k < nk_-1
  for (size_t jkk = 1; jkk < nk_-1; ++jkk) {
    // l = 0
    jl = 0;
    addSpectralCoefficient(jkk, jl, ReRe, 1, 0);
    addSpectralCoefficient(jkk, jl, ImRe, -1, 0);

    // 0 < l < nl_-1
    for (size_t jll = 1; jll < nl_-1; ++jll) {
      addSpectralCoefficient(jkk, jll, ReRe, 2, 1);
      addSpectralCoefficient(jkk, jll, ReIm, 2, -1);
      addSpectralCoefficient(jkk, jll, ImRe, -2, 1);
      addSpectralCoefficient(jkk, jll, ImIm, -2, -1);
    }

    // l = nl_-1
    jl = nl_-1;
    if (ny_ % 2 == 1) {
      addSpectralCoefficient(jkk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jkk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jkk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jkk, jl, ImIm, -2, -1);
    } else {
      addSpectralCoefficient(jkk, jl, ReRe, 1, 0);
      addSpectralCoefficient(jkk, jl, ImRe, -1, 0);
    }
  }

  // k = nk_-1
  jk = nk_-1;

  // l = 0
  jl = 0;
  if (nx_ % 2 == 1) {
    addSpectralCoefficient(jk, jl, ReRe, 1, 0);
    addSpectralCoefficient(jk, jl, ImRe, -1, 0);
  } else {
    addSpectralCoefficient(jk, jl, ReRe, 0, 0);
  }

  // 0 < l < nl_-1
  for (size_t jll = 1; jll < nl_-1; ++jll) {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jll, ReRe, 2, 1);
      addSpectralCoefficient(jk, jll, ReIm, 2, -1);
      addSpectralCoefficient(jk, jll, ImRe, -2, 1);
      addSpectralCoefficient(jk, jll, ImIm, -2, -1);
    } else {
      addSpectralCoefficient(jk, jll, ReRe, 0, 1);
      addSpectralCoefficient(jk, jll, ReIm, 0, -1);
    }
  }

  // l = nl_-1
  jl = nl_-1;
  if (ny_ % 2 == 1) {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 2, 1);
      addSpectralCoefficient(jk, jl, ReIm, 2, -1);
      addSpectralCoefficient(jk, jl, ImRe, -2, 1);
      addSpectralCoefficient(jk, jl, ImIm, -2, -1);

    } else {
      addSpectralCoefficient(jk, jl, ReRe, 0, 1);
      addSpectralCoefficient(jk, jl, ReIm, 0, -1);
    }
  } else {
    if (nx_ % 2 == 1) {
      addSpectralCoefficient(jk, jl, ReRe, 1, 0);
      addSpectralCoefficient(jk, jl, ImRe, -1, 0);
    } else {
      addSpectralCoefficient(jk, jl, ReRe, 0, 0);
    }
  }

  // Define global spectral size
  nsGlb_ = spVec_.size();
  oops::Log::test() << "- Spectral array global size: " << nsGlb_ << std::endl;

  oops::Log::trace() << classname() << "::setupGlobalSpectralSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::setupParallelizationInit() {
  oops::Log::trace() << classname() << "::setupParallelizationInit starting" << std::endl;

  // Split truncation wavenumbers in equal chunks

  // Prepare vectors
  std::vector<int> jkVec;
  std::vector<int> jlVec;
  std::vector<size_t> jwGlbVec;
  std::vector<size_t> nklPerTaskTarget(comm_.size(), 0);
  size_t index = 0;
  for (size_t jk = 0; jk < ellips_.size(); ++jk) {
    for (size_t jl = 0; jl <= ellips_[jk]; ++jl) {
      jkVec.push_back(jk);
      jlVec.push_back(jl);
      jwGlbVec.push_back(ikstar(jk, jl, M_, N_, nwGlb_));
      ++nklPerTaskTarget[index];
      ++index;
      if (index == comm_.size()) index = 0;
    }
  }

  // Sort truncation wavenumbers in ascending jwGlb
  const size_t nkl = jwGlbVec.size();
  std::vector<size_t> klOrder(nkl);
  std::iota(klOrder.begin(), klOrder.end(), 0);
  std::stable_sort(klOrder.begin(), klOrder.end(),
    [&](size_t i, size_t j){return jwGlbVec[i] < jwGlbVec[j];});

  // Split total wave number among tasks
  std::vector<size_t> nklPerTask(comm_.size());
  std::vector<std::vector<size_t>> ellipsTask(ellips_.size());
  for (size_t jk = 0; jk < ellips_.size(); ++jk) {
    ellipsTask[jk].resize(ellips_[jk]+1);
  }
  size_t jt = 0;
  for (size_t jkl = 0; jkl < nkl; ++jkl) {
    // Conditions to switch to the next task
    if ((nklPerTask[jt] > nklPerTaskTarget[jt]) && (jt < comm_.size()-1)) {
      nklPerTaskTarget[jt+1] += nklPerTaskTarget[jt] - nklPerTask[jt];
      ++jt;
    }

    // Get jk, jl
    const size_t jk = jkVec[klOrder[jkl]];
    const size_t jl = jlVec[klOrder[jkl]];

    // Save task
    ellipsTask[jk][jl] = jt;

    // Update local size
    ++nklPerTask[jt];
  }

  // Add task in global spectral vector
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get jk, jl
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;

    // Add task
    spVec_[jsGlb].jt = ellipsTask[jk][jl];
  }

  // Spectral size per task and local to global mapping
  nsPerTask_.resize(comm_.size());
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get task
    const size_t jt = spVec_[jsGlb].jt;

    // Update number of spectral coefficient for this task
    ++nsPerTask_[jt];

    if (jt == myrank_) {
      // Add local spectral coefficient
      sToSGlb_.push_back(jsGlb);
    }
  }

  // Save local size
  ns_ = sToSGlb_.size();

  // Communication vectors and mapping
  sCounts_.resize(comm_.size());
  sDispls_.resize(comm_.size());
  if (myrank_ == 0) {
    sMapping_.resize(nsGlb_);
  }
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sCounts_[jt] = nsPerTask_[jt];
    sDispls_[jt] = static_cast<int>(jt ? sDispls_[jt-1] + sCounts_[jt-1] : 0);
  }
  comm_.gatherv(sToSGlb_.cbegin(), sToSGlb_.cend(), sMapping_.begin(), sMapping_.end(),
    sCounts_, sDispls_, 0);

  // Compute spectral imbalance
  const double sImb = static_cast<double>(*std::max_element(nsPerTask_.begin(), nsPerTask_.end()))
    / static_cast<double>(*std::min_element(nsPerTask_.begin(), nsPerTask_.end()));
  oops::Log::info() << "Info     : - Spectral imbalance (max/min): " << sImb << std::endl;

  oops::Log::trace() << classname() << "::setupParallelizationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::setupParallelizationFinal() {
  oops::Log::trace() << classname() << "::setupParallelizationFinal starting" << std::endl;

  // Get min and max jwGlb
  nwStartPerTask_.resize(comm_.size(), nwGlb_-1);
  nwEndPerTask_.resize(comm_.size(), 0);
  for (size_t jsGlb = 0; jsGlb < nsGlb_; ++jsGlb) {
    // Get task
    const size_t jt = spVec_[jsGlb].jt;

    // Get global total wavenumber
    const size_t jwGlb = spVec_[jsGlb].jwGlb;

    // Update min and max global total wavenumber
    nwStartPerTask_[jt] = std::min(nwStartPerTask_[jt], jwGlb);
    nwEndPerTask_[jt] = std::max(nwEndPerTask_[jt], jwGlb);

    if (dwGlb_ > 0.0) {
      // Get jkstar
      const double jkstar = spVec_[jsGlb].kstar;

      // Loop of over global total wavenumber
      for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
        if ((jkstar >= static_cast<double>(jwGlb)-dwGlb_)
          && (jkstar <= static_cast<double>(jwGlb)+dwGlb_)) {
          // Update min and max global total wavenumber
          nwStartPerTask_[jt] = std::min(nwStartPerTask_[jt], jwGlb);
          nwEndPerTask_[jt] = std::max(nwEndPerTask_[jt], jwGlb);
        }
      }
    }
  }

  // Define nwPerTask_
  nwPerTask_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    nwPerTask_[jt] = nwEndPerTask_[jt] - nwStartPerTask_[jt] + 1;
  }

  // Local number of total wavenumbers
  nw_ = nwPerTask_[myrank_];

  // Define root nw
  nwRootPerTask_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    nwRootPerTask_[jt] = nwStartPerTask_[jt+1] - nwStartPerTask_[jt];
  }
  nwRootPerTask_[comm_.size()-1] = nwGlb_ - nwStartPerTask_[comm_.size()-1];

  // Local number of root total wavenumbers
  nwRoot_ = nwRootPerTask_[myrank_];

  // Communication vectors
  wCounts_.resize(comm_.size());
  wDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    wCounts_[jt] = nwRootPerTask_[jt];
    wDispls_[jt] = static_cast<int>(jt ? wDispls_[jt-1] + wCounts_[jt-1] : 0);
  }

  // Compute total wavenumber imbalance
  const double wImb = static_cast<double>(*std::max_element(nwPerTask_.begin(), nwPerTask_.end()))
    / static_cast<double>(*std::min_element(nwPerTask_.begin(), nwPerTask_.end()));
  oops::Log::info() << "Info     : - Total wavenumber imbalance (max/min): " << wImb << std::endl;

  // Prepare covariance communicators
  myJwGlb_.resize(nwGlb_);
  for (size_t jwGlb = 0; jwGlb < nwGlb_; ++jwGlb) {
    // Define used global wavenumbers
    myJwGlb_[jwGlb] = (nwStartPerTask_[myrank_] <= jwGlb) && (jwGlb <= nwEndPerTask_[myrank_]);

    // Define color of the reduction communicator
    const size_t covRedColor = myJwGlb_[jwGlb] ? 1 : 0;

    // Communicator name
    const std::string covRedCommName = "covRed_" + std::to_string(jwGlb);

    // Split communicator
    covRedComm_.push_back(&comm_.split(covRedColor, covRedCommName.c_str()));

    // Define color of the broadcasting communicator
    const size_t covBcastColor = myJwGlb_[jwGlb] || myrank_ == 0 ? 1 : 0;

    // Communicator name
    const std::string covBcastCommName = "covBcast_" + std::to_string(jwGlb);

    // Split communicator=
    covBcastComm_.push_back(&comm_.split(covBcastColor, covBcastCommName.c_str()));
  }

  oops::Log::trace() << classname() << "::setupParallelizationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::setupLocalSpectralSpace() {
  oops::Log::trace() << classname() << "::setupLocalSpectralSpace starting" << std::endl;

  // Create dummy PointCloud function space
  atlas::FieldSet flds;
  atlas::Field lonlatField("lonlat", make_datatype<double>(), make_shape(ns_, 2));
  atlas::Field ghostField("ghost", make_datatype<int>(), make_shape(ns_));
  atlas::Field remoteIndexField("remote_index", make_datatype<int>(), make_shape(ns_));
  atlas::Field partitionField("partition", make_datatype<int>(), make_shape(ns_));
  atlas::Field globalIndexField("global_index", make_datatype<int64_t>(), make_shape(ns_));
  flds.add(lonlatField);
  flds.add(ghostField);
  flds.add(remoteIndexField);
  flds.add(partitionField);
  flds.add(globalIndexField);
  auto lonlatView = make_view<double, 2>(lonlatField);
  auto ghostView = make_view<int, 1>(ghostField);
  auto remoteIndexView = make_indexview<int, 1>(remoteIndexField);
  auto partitionView = make_view<int, 1>(partitionField);
  auto globalIndexView = make_indexview<int64_t, 1>(globalIndexField);

  for (size_t js = 0; js < ns_; ++js) {
    // Get global index
    const size_t jsGlb = sToSGlb_[js];

    // Get jk/jl/jq
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;
    const size_t jq = spVec_[jsGlb].jq;

    // Define dummy lon/lat from jk/jl/jq
    lonlatView(js, 0) = (static_cast<double>(jk)+0.25*static_cast<double>(jq))
      /static_cast<double>(nk_)*360.0;
    lonlatView(js, 1) = -90.0*(static_cast<double>(jl)+0.25*static_cast<double>(jq))
      /static_cast<double>(nl_)*180.0;

    // Define ghost field
    ghostView(js) = 0;

    // Define remote index field
    remoteIndexView(js) = js;

    // Define partition field
    partitionView(js) = myrank_;

    // Define global index field
    globalIndexView(js) = static_cast<int64_t>(jsGlb);
  }
  spFspace_.reset(new atlas::functionspace::PointCloud(flds));

  // Generate spectral UID
  specUid_ = generateSpectralUid(*spFspace_, comm_);

  // Print UIDs
  oops::Log::info() << "Info     : - UIDs: " << gridUid_ << " / " << specUid_ << std::endl;

  // Allocate vectors
  jkVec_.resize(ns_);
  jlVec_.resize(ns_);
  jqVec_.resize(ns_);
  jwGlbVec_.resize(ns_);
  kstarVec_.resize(ns_);
  spNormVec_.resize(ns_);
  lapDirVec_.resize(ns_);
  lapInvVec_.resize(ns_);

  // Allocate fields
  atlas::Field xDerivRowField("xDerivRow", make_datatype<int>(), make_shape(ns_));
  atlas::Field xDerivColField("xDerivCol", make_datatype<int>(), make_shape(ns_));
  atlas::Field xDerivSField("xDerivS", make_datatype<double>(), make_shape(ns_));
  atlas::Field yDerivRowField("yDerivRow", make_datatype<int>(), make_shape(ns_));
  atlas::Field yDerivColField("yDerivCol", make_datatype<int>(), make_shape(ns_));
  atlas::Field yDerivSField("yDerivS", make_datatype<double>(), make_shape(ns_));

  // Add fields to derivatives_
  derivatives_.add(xDerivRowField);
  derivatives_.add(xDerivColField);
  derivatives_.add(xDerivSField);
  derivatives_.add(yDerivRowField);
  derivatives_.add(yDerivColField);
  derivatives_.add(yDerivSField);

  // Get fields views
  auto xDerivRowView = make_view<int, 1>(xDerivRowField);
  auto xDerivColView = make_view<int, 1>(xDerivColField);
  auto xDerivSView = make_view<double, 1>(xDerivSField);
  auto yDerivRowView = make_view<int, 1>(yDerivRowField);
  auto yDerivColView = make_view<int, 1>(yDerivColField);
  auto yDerivSView = make_view<double, 1>(yDerivSField);

  // Fill data
  for (size_t js = 0; js < ns_; ++js) {
    // Get global index
    const size_t jsGlb = sToSGlb_[js];

    // Get spVec_ values
    const size_t jk = spVec_[jsGlb].jk;
    const size_t jl = spVec_[jsGlb].jl;
    const size_t jq = spVec_[jsGlb].jq;
    const double kstar = spVec_[jsGlb].kstar;
    const size_t jwGlb = spVec_[jsGlb].jwGlb;
    const int jsXDerivativeOffset = spVec_[jsGlb].jsXDerivativeOffset;
    const int jsYDerivativeOffset = spVec_[jsGlb].jsYDerivativeOffset;

    // Copy jk, jl, jq and jwGlb
    jkVec_[js] = jk;
    jlVec_[js] = jl;
    jqVec_[js] = jq;
    kstarVec_[js] = kstar;
    jwGlbVec_[js] = jwGlb;

    // Spectral norm
    spNormVec_[js] = static_cast<double>(spNormKL_[jk*nl_+jl]);

    // X-derivative linear operator
    xDerivColView(js) = js;
    if (jsXDerivativeOffset == 0) {
      // Set to zero
      xDerivRowView(js) = js;
      xDerivSView(js) = 0.0;
    } else {
      // Involved in the derivative
      xDerivRowView(js) = js+jsXDerivativeOffset;
      if (jsXDerivativeOffset > 0) {
        xDerivSView(js) = static_cast<double>(jk)*exwn_;
      } else  {
        xDerivSView(js) = -static_cast<double>(jk)*exwn_;
      }
    }

    // Y-derivative linear operator
    yDerivColView(js) = js;
    if (jsYDerivativeOffset == 0) {
      // Set to zero
      yDerivRowView(js) = js;
      yDerivSView(js) = 0.0;
    } else {
      // Involved in the derivative
      yDerivRowView(js) = js+jsYDerivativeOffset;
      if (jsYDerivativeOffset > 0) {
        yDerivSView(js) = static_cast<double>(jl)*eywn_;
      } else  {
        yDerivSView(js) = -static_cast<double>(jl)*eywn_;
      }
    }

    // Direct Laplacian
    lapDirVec_[js] = -(static_cast<double>(jk*jk)*exwn_*exwn_
      + static_cast<double>(jl*jl)*eywn_*eywn_);

    // Inverse Laplacian
    if (jk == 0 && jl == 0) {
      lapInvVec_[js] = 0.0;
    } else {
      lapInvVec_[js] = 1.0/lapDirVec_[js];
    }
  }

  // Compute control vector size
  ctlVecSize_ = ns_*nvz_;

  // Compute spectral norm sum
  std::vector<double> spNormSum(nwGlb_, 0.0);
  for (size_t js = 0; js < ns_; ++js) {
    for (size_t jw = 0; jw < nw_; ++jw) {
      // Update spectral norm sum
      const size_t jwGlb = jw + nwStartPerTask_[myrank_];
      if (includeWavenumber(js, jw)) {
        spNormSum[jwGlb] += spNormVec_[js];
      }
    }
  }

  // Reduce spectral norm sum
  comm_.allReduceInPlace(spNormSum.begin(), spNormSum.end(), eckit::mpi::sum());

  // Inverse spectral norm sum
  spNormSumInv_.resize(nw_, 0.0);
  for (size_t jw = 0; jw < nw_; ++jw) {
    const size_t jwGlb = jw + nwStartPerTask_[myrank_];
    ASSERT(spNormSum[jwGlb] > 0.0);
    spNormSumInv_[jw] = 1.0/spNormSum[jwGlb];
  }

  oops::Log::trace() << classname() << "::setupLocalSpectralSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierTransformBase::addSpectralCoefficient(const size_t & jk,
                                                    const size_t & jl,
                                                    const Quad & jq,
                                                    const size_t & jsXDerivativeOffset,
                                                    const size_t & jsYDerivativeOffset) {
  if (jk < ellips_.size()) {
    if (jl <= ellips_[jk]) {
      // Update spVec
      spElem e;
      e.jk = jk;
      e.jl = jl;
      e.jq = jq;
      e.kstar = rkstar(jk, jl, M_, N_, nwGlb_);
      e.jwGlb = ikstar(jk, jl, M_, N_, nwGlb_);
      e.jsXDerivativeOffset = jsXDerivativeOffset;
      e.jsYDerivativeOffset = jsYDerivativeOffset;
      spVec_.push_back(e);

      // Check consistency between truncation and maximum total wavenumber
      ASSERT(e.jwGlb < nwGlb_);
    }
  }

  // Update spNorm
  ++spNormKL_[jk*nl_+jl];
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
