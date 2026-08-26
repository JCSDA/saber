/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierSpectralConverter.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "oops/util/FloatCompare.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierSpectralConverter>
  makerBifourierSpectralConverter_("BifourierSpectralConverter");

// -----------------------------------------------------------------------------

BifourierSpectralConverter::BifourierSpectralConverter(const oops::GeometryData & outerGeometryData,
                                                       const oops::Variables & outerVars,
                                                       const eckit::Configuration & covarConfig,
                                                       const Parameters_ & params,
                                                       const oops::FieldSet3D & xb,
                                                       const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    comm_(outerGeometryData.comm()),
    myrank_(comm_.rank()),
    innerVars_(outerVars)
{
  oops::Log::trace() << classname() << "::BifourierSpectralConverter starting" << std::endl;

  // Retrieve outer spectral transform
  outerTrans_ = transStore_.retrieveTransform(outerGeometryData, outerVars);

  // Create inner grid-point FunctionSpace
  atlas::functionspace::StructuredColumns innerGpFs;
  if (params.fspaceFromBkgVar.value()) {
    // Use the function space of a field of the background
    innerGpFs = xb[*params.fspaceFromBkgVar.value()].functionspace();
  } else {
    // Get outer geometry configuration
    const atlas::functionspace::StructuredColumns outerFs(
      outerTrans_->geometryData().functionSpace());
    const atlas::StructuredGrid & outerGrid = outerFs.grid();
    const atlas::util::Config outerGridConfig = outerGrid.spec();

    // Get xSpace and ySpace properties
    atlas::util::Config xSpaceConfig = outerGridConfig.getSubConfiguration("xspace");
    atlas::util::Config ySpaceConfig = outerGridConfig.getSubConfiguration("yspace");
    const int outerNx = xSpaceConfig.getInt("N");
    const int outerNy = ySpaceConfig.getInt("N");
    const double outerStartX = xSpaceConfig.getDouble("start");
    const double outerStartY = ySpaceConfig.getDouble("start");
    const double outerEndX = xSpaceConfig.getDouble("end");
    const double outerEndY = ySpaceConfig.getDouble("end");
    const double outerDx = (outerEndX-outerStartX)/static_cast<double>(outerNx-1);
    const double outerDy = (outerEndY-outerStartY)/static_cast<double>(outerNy-1);

    // Get domain size, assuming a periodic domain
    const double Lx = static_cast<double>(outerNx)*outerDx;
    const double Ly = static_cast<double>(outerNy)*outerDy;

    // Update xSpace and ySpace properties
    ASSERT(params.nx.value() && params.ny.value());
    const int innerNx = *params.nx.value();
    const int innerNy = *params.ny.value();
    const double innerDx = Lx/static_cast<double>(innerNx);
    const double innerDy = Ly/static_cast<double>(innerNy);
    const double innerStartX = outerStartX;
    const double innerStartY = outerStartY;
    const double innerEndX = innerStartX+innerDx*static_cast<double>(innerNx-1);
    const double innerEndY = innerStartY+innerDy*static_cast<double>(innerNy-1);

    // Check consistency
    const double ratioNx = static_cast<double>(outerNx)/static_cast<double>(innerNx);
    const double ratioNy = static_cast<double>(outerNy)/static_cast<double>(innerNy);
    ASSERT(oops::is_close_relative(ratioNx, ratioNy, 1.0e-12));

    // Update xSpace and ySpace configurations
    xSpaceConfig.set("N", innerNx);
    ySpaceConfig.set("N", innerNy);
    xSpaceConfig.set("start", innerStartX);
    ySpaceConfig.set("start", innerStartY);
    xSpaceConfig.set("end", innerEndX);
    ySpaceConfig.set("end", innerEndY);

    // Get and update domain configuration
    atlas::util::Config domainConfig = outerGridConfig.getSubConfiguration("domain");
    domainConfig.set("xmin", innerStartX);
    domainConfig.set("ymin", innerStartY);
    domainConfig.set("xmax", innerEndX);
    domainConfig.set("ymax", innerEndY);

    // Get projection (same for outer and inner geometries)
    const atlas::Projection projection(outerGridConfig.getSubConfiguration("projection"));

    // Create new xSpace and ySpace
    const atlas::StructuredGrid::XSpace xspace(xSpaceConfig);
    const atlas::StructuredGrid::YSpace yspace(ySpaceConfig);

    // Create new domain
    const atlas::Domain domain(domainConfig);

    // Create inner grid
    atlas::StructuredGrid innerGrid(xspace, yspace, projection, domain);

    // Set up ATLAS MPI
    eckit::mpi::setCommDefault(comm_.name().c_str());

    // Create partitioner
    atlas::grid::Partitioner partitioner = atlas::grid::Partitioner(params.partitioner.value());

    // Create functionspace from partitioner
    innerGpFs = atlas::functionspace::StructuredColumns(innerGrid, partitioner);

    // Reset communicator
    eckit::mpi::setCommDefault(comm_.name().c_str());

    // Bugfix for regional grids
    // It seems that the content of lonlat for a regional function space is actually the xy
    // coordinates. The routine to compute distances on the sphere was complaining about
    // impossible lon/lat values...
    auto lonlat = atlas::array::make_view<double, 2>(innerGpFs.lonlat());
    double lonlatPoint[] = {0, 0};
    const auto view_i = atlas::array::make_indexview<int, 1>(innerGpFs.index_i());
    const auto view_j = atlas::array::make_indexview<int, 1>(innerGpFs.index_j());
    for (int jj = 0; jj < innerGpFs.size(); ++jj) {
      innerGrid.lonlat(view_i(jj), view_j(jj), lonlatPoint);
      lonlat(jj, 0) = lonlatPoint[0];
      lonlat(jj, 1) = lonlatPoint[1];
    }
  }

  // Copy 1D fields (profiles), no interpolation for gridded fields
  atlas::FieldSet fields;
  for (const auto & field : outerGeometryData.fieldSet()) {
    if (field.rank() == 1) {
      fields.add(field);
    }
  }

  // Inner geometry data
  innerGpGeometryData_ = std::make_unique<oops::GeometryData>(innerGpFs,
    fields, outerGeometryData.levelsAreTopDown(), comm_, false);

  // Create inner spectral transform
  innerTrans_ = transStore_.setupTransform(*innerGpGeometryData_, innerVars_,
    params.transform.value());

  // Create inner spectral GeometryData
  innerGeometryData_ = std::make_unique<oops::GeometryData>(innerTrans_->spFspace(),
    fields, outerGeometryData.levelsAreTopDown(), comm_, false);

  // Check domain size
  ASSERT(oops::is_close_relative(static_cast<double>(innerTrans_->nx())*innerTrans_->dx(),
    static_cast<double>(outerTrans_->nx())*outerTrans_->dx(), 1.0e-12));
  ASSERT(oops::is_close_relative(static_cast<double>(innerTrans_->ny())*innerTrans_->dy(),
    static_cast<double>(outerTrans_->ny())*outerTrans_->dy(), 1.0e-12));

  // Prepare spectral converter mapping
  std::vector<int> outerToInnerJsGlb(outerTrans_->nsGlb(), -1);
  if (outerTrans_->nsGlb() >= innerTrans_->nsGlb()) {
    // Zero-padding case
    size_t outerJsGlb = 0;
    for (size_t innerJsGlb = 0; innerJsGlb < innerTrans_->nsGlb(); ++innerJsGlb) {
      // Increment outer jsGlb
      while ((innerTrans_->sGlbToK(innerJsGlb) != outerTrans_->sGlbToK(outerJsGlb)) &&
        (innerTrans_->sGlbToL(innerJsGlb) != outerTrans_->sGlbToL(outerJsGlb))) {
        ++outerJsGlb;
      }

      // Check jq
      if (innerTrans_->sGlbToQ(innerJsGlb) == outerTrans_->sGlbToQ(outerJsGlb)) {
        ASSERT(outerToInnerJsGlb[outerJsGlb] == -1);
        outerToInnerJsGlb[outerJsGlb] = innerJsGlb;
      }

      // Increment outer jsGlb
      ++outerJsGlb;
    }
  } else {
    // Truncation case
    size_t innerJsGlb = 0;
    for (size_t outerJsGlb = 0; outerJsGlb < outerTrans_->nsGlb(); ++outerJsGlb) {
      // Increment inner jsGlb
      while ((innerTrans_->sGlbToK(innerJsGlb) != outerTrans_->sGlbToK(outerJsGlb)) &&
        (innerTrans_->sGlbToL(innerJsGlb) != outerTrans_->sGlbToL(outerJsGlb))) {
        ++innerJsGlb;
      }

      // Check jq
      if (innerTrans_->sGlbToQ(innerJsGlb) == outerTrans_->sGlbToQ(outerJsGlb)) {
        ASSERT(outerToInnerJsGlb[outerJsGlb] == -1);
        outerToInnerJsGlb[outerJsGlb] = innerJsGlb;
      }

      // Increment inner jsGlb
      ++innerJsGlb;
    }
  }

  // Prepare spectral converter communications
  std::vector<int> sendOuterIndexGlb;
  sendCounts_.resize(comm_.size());
  recvCounts_.resize(comm_.size());
  std::fill(sendCounts_.begin(), sendCounts_.end(), 0);
  std::fill(recvCounts_.begin(), recvCounts_.end(), 0);
  for (size_t outerJsGlb = 0; outerJsGlb < outerTrans_->nsGlb(); ++outerJsGlb) {
    // Get ordered index
    const size_t orderedOuterJsGlb = outerTrans_->sMapping()[outerJsGlb];

    if (outerToInnerJsGlb[orderedOuterJsGlb] >= 0) {
      // Get inner jsGlb
      const size_t innerJsGlb = outerToInnerJsGlb[orderedOuterJsGlb];

      // Get inner and outer tasks
      const size_t innerTask = innerTrans_->sGlbToTask(innerJsGlb);
      const size_t outerTask = outerTrans_->sGlbToTask(orderedOuterJsGlb);

      if (innerTask == myrank_) {
        // Save coefficients to send
        const size_t innerJs = innerTrans_->sGlbToS(innerJsGlb);
        sendIndex_.push_back(innerJs);
        sendOuterIndexGlb.push_back(orderedOuterJsGlb);

        // Update send counts
        ++sendCounts_[outerTask];
      }

      if (outerTask == myrank_) {
        // Update receive counts
        ++recvCounts_[innerTask];
      }
    }
  }

  // Sizes
  sendSize_ = 0;
  for (const auto & n : sendCounts_) sendSize_ += n;
  recvSize_ = 0;
  for (const auto & n : recvCounts_) recvSize_ += n;

  // Displacements
  sendDispls_.resize(comm_.size());
  recvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendDispls_[jt] = static_cast<int>(jt ? sendDispls_[jt-1] + sendCounts_[jt-1] : 0);
    recvDispls_[jt] = static_cast<int>(jt ? recvDispls_[jt-1] + recvCounts_[jt-1] : 0);
  }

  // Communicate receive global index
  std::vector<int> recvOuterIndexGlb(recvSize_);
  comm_.allToAllv(sendOuterIndexGlb.data(), sendCounts_.data(), sendDispls_.data(),
    recvOuterIndexGlb.data(), recvCounts_.data(), recvDispls_.data());

  // Retrieve local receive index
  recvIndex_.resize(recvSize_);
  for (size_t jj = 0; jj < recvSize_; ++jj) {
    recvIndex_[jj] = outerTrans_->sGlbToS(recvOuterIndexGlb[jj]);
  }

  // Multiply by nvz
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendCounts_[jt] *= innerTrans_->nvz();
    recvCounts_[jt] *= outerTrans_->nvz();
    sendDispls_[jt] *= innerTrans_->nvz();
    recvDispls_[jt] *= outerTrans_->nvz();
  }

  oops::Log::trace() << classname() << "::BifourierSpectralConverter done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralConverter::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Create buffers
  std::vector<double> sendVec(sendSize_*innerTrans_->nvz(), 0.0);
  std::vector<double> recvVec(recvSize_*outerTrans_->nvz());

  // Reserialize from spectral FieldSet
  size_t zOffset = 0;
  for (const auto & var : outerVars_) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = fset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(innerTrans_->ns()));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jj = 0; jj < sendSize_; ++jj) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = sendIndex_[jj];

        // Communication vector index
        const size_t jjv = jj*innerTrans_->nvz()+jvz;

        // Copy data
        sendVec[jjv] = spView(js, jz);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  comm_.allToAllv(sendVec.data(), sendCounts_.data(), sendDispls_.data(),
    recvVec.data(), recvCounts_.data(), recvDispls_.data());

  // Remove outer variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), outerVars_.variables());

  // Deserialize into grid-point FieldSet
  zOffset = 0;
  for (const auto & var : outerVars_) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Create field
    atlas::Field spField = outerGeometryData_.functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    fset.add(spField);

    // Get field view
    auto spView = make_view<double, 2>(spField);

    // Set to zero
    spView.assign(0.0);

    for (size_t jj = 0; jj < recvSize_; ++jj) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = recvIndex_[jj];

        // Communication vector index
        const size_t jjv = jj*outerTrans_->nvz()+jvz;

        // Copy data
        spView(js, jz) = recvVec[jjv];
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralConverter::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Create buffers
  std::vector<double> recvVec(recvSize_*outerTrans_->nvz());
  std::vector<double> sendVec(sendSize_*innerTrans_->nvz());

  // Reserialize from spectral FieldSet
  size_t zOffset = 0;
  for (const auto & var : outerVars_) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Check field
    const auto spField = fset[var.name()];
    ASSERT(spField.shape(0) == static_cast<int>(outerTrans_->ns()));
    ASSERT(spField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto spView = make_view<double, 2>(spField);

    for (size_t jj = 0; jj < recvSize_; ++jj) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = recvIndex_[jj];

        // Communication vector index
        const size_t jjv = jj*outerTrans_->nvz()+jvz;

        // Copy data
        recvVec[jjv] = spView(js, jz);
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  // Communication
  comm_.allToAllv(recvVec.data(), recvCounts_.data(), recvDispls_.data(),
    sendVec.data(), sendCounts_.data(), sendDispls_.data());

  // Remove inner variables
  util::removeFieldsFromFieldSet(fset.fieldSet(), innerVars_.variables());

  // Deserialize into grid-point FieldSet
  zOffset = 0;
  for (const auto & var : outerVars_) {
    // Get number of vertical levels
    const size_t nz = var.getLevels();

    // Create field
    atlas::Field spField = innerGeometryData_->functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    fset.add(spField);

    // Get field view
    auto spView = make_view<double, 2>(spField);

    // Set to zero
    spView.assign(0.0);

    for (size_t jj = 0; jj < sendSize_; ++jj) {
      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Spectral index
        const size_t js = sendIndex_[jj];

        // Communication vector index
        const size_t jjv = jj*innerTrans_->nvz()+jvz;

        // Copy data
        spView(js, jz) = sendVec[jjv];
      }
    }

    // Update total number of levels
    zOffset += nz;
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralConverter::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  // Check sizes
  ASSERT(innerTrans_->nx() <= outerTrans_->nx());
  ASSERT(innerTrans_->ny() <= outerTrans_->ny());

  // Adjoint multiply
  multiplyAD(fset);

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierSpectralConverter::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
