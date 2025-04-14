/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BiperiodizationImpl.h"

#include <algorithm>

#include "atlas/grid.h"

#include "oops/util/FloatCompare.h"

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

BiperiodizationImpl::BiperiodizationImpl(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & vars,
                                         const Parameters_ & params)
  : outerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    myrank_(comm_.rank()),
    vars_(vars),
    params_(params)
{
  oops::Log::trace() << classname() << "::BiperiodizationImpl starting" << std::endl;

  // Check function space type
  ASSERT(outerGeometryData.functionSpace().type() == "StructuredColumns");

  // Number of levels for all variables
  nvz_ = 0;
  for (const auto & var : vars_) {
    nvz_ += var.getLevels();
  }

  // Get outer grid xspace and yspace
  const atlas::functionspace::StructuredColumns outerFs(outerGeometryData.functionSpace());
  const auto outerGrid = outerFs.grid();
  atlas::util::Config xspec = outerGrid.xspace().spec();
  atlas::util::Config yspec = outerGrid.yspace().spec();
  const size_t outerNx = xspec.getInt("N");
  const double startx = xspec.getDouble("start");
  const double outerEndx = xspec.getDouble("end");
  const double dx = (outerEndx-startx)/static_cast<double>(outerNx-1);
  const size_t outerNy = yspec.getInt("N");
  const double starty = yspec.getDouble("start");
  const double outerEndy = yspec.getDouble("end");
  const double dy = (outerEndy-starty)/static_cast<double>(outerNy-1);
  oops::Log::info() << "Info     : Outer grid: " << std::endl;
  oops::Log::test() << "- xspace: " << xspec << std::endl;
  oops::Log::test() << "- yspace: " << yspec << std::endl;
  oops::Log::test() << "Outer grid size: " << outerGrid.size() << std::endl;

  // Get inner and outer extension zone
  const size_t innerExtNx = params_.innerExtNx.value();
  const size_t innerExtNy = params_.innerExtNy.value();
  const size_t outerExtNx = params_.outerExtNx.value();
  const size_t outerExtNy = params_.outerExtNy.value();

  // Define physical grid size
  const size_t physicalNx = outerNx - outerExtNx;
  const size_t physicalNy = outerNy - outerExtNy;

  if (innerExtNx == outerExtNx && innerExtNy == outerExtNy) {
    // Copy grid
    innerGrid_ = outerGrid;
    oops::Log::info() << "Inner grid = outer grid" << std::endl;
  } else {
    // Define inner grid
    const size_t innerNx = physicalNx + innerExtNx;
    const size_t innerNy = physicalNy + innerExtNy;
    const double innerEndx = startx+static_cast<double>(innerNx-1)*dx;
    const double innerEndy = starty+static_cast<double>(innerNy-1)*dy;
    xspec.set("end", innerEndx);
    yspec.set("end", innerEndy);
    xspec.set("N", innerNx);
    yspec.set("N", innerNy);
    oops::Log::info() << "Info     : Inner grid: " << std::endl;
    oops::Log::test() << "- xspace: " << xspec << std::endl;
    oops::Log::test() << "- yspace: " << yspec << std::endl;

    // Generate inner StructuredGrid
    atlas::grid::detail::grid::Structured::XSpace innerXSpace(xspec);
    atlas::grid::detail::grid::Structured::YSpace innerYSpace(yspec);
    innerGrid_ = atlas::StructuredGrid(innerXSpace, innerYSpace, outerFs.grid().projection());
    oops::Log::test() << "Inner grid size: " << innerGrid_.size() << std::endl;
  }

  // Allocate inner partition
  innerPartition_.resize(innerGrid_.size());

  // Get outer Partitioner name
  const std::string outerPartitionerName = outerFs.distribution();
  oops::Log::info() << "Outer partitioner: " << outerPartitionerName << std::endl;

  if (outerPartitionerName == "custom") {
    // Get inner Partitioner
    const std::string innerPartitionerName = params_.innerPartitioner.value();
    oops::Log::info() << "Inner partitioner: " << innerPartitionerName << std::endl;
    const auto innerPartitioner =  atlas::grid::Partitioner(innerPartitionerName);

    // Generate inner partition
    innerPartitioner.partition(innerGrid_, innerPartition_.data());
  } else {
    // Get outer Partitioner
    const auto outerPartitioner =  atlas::grid::Partitioner(outerPartitionerName);
    oops::Log::info() << "Inner partitioner = outer partitioner ("
       << outerFs.distribution() << ")" << std::endl;

    // Generate inner partition
    outerPartitioner.partition(innerGrid_, innerPartition_.data());
  }

  if (innerExtNx == outerExtNx && innerExtNy == outerExtNy && outerPartitionerName != "custom") {
    // Copy function space
    innerFs_ = outerFs;
  } else {
    // Generate inner Distribution
    atlas::grid::Distribution innerDistribution(comm_.size(), innerGrid_.size(),
    innerPartition_.data());

    // Generate inner FunctionSpace
    innerFs_ = atlas::functionspace::StructuredColumns(innerGrid_, innerDistribution);
  }

  // Prepare mixing mask components
  const size_t nmix = params_.nmix.value() != boost::none ? *params_.nmix.value() :
    std::max(outerExtNx, outerExtNy);
  const double Lmix = params_.Lmix.value();
  std::vector<double> mixingX(outerExtNx);
  for (size_t jx = 0; jx < outerExtNx; ++jx) {
    const double ux = static_cast<double>(jx+1)/static_cast<double>(std::min(outerExtNx, nmix)+1);
    if (ux < 1.0) {
      mixingX[jx] = 0.5*(1.0+std::erf(Lmix*(1.0-2.0*ux)/std::sqrt(4.0*ux*(1.0-ux))));
    } else {
      mixingX[jx] = 0.0;
    }
  }
  std::vector<double> mixingY(outerExtNy);
  for (size_t jy = 0; jy < outerExtNy; ++jy) {
    const double uy = static_cast<double>(jy+1)/static_cast<double>(std::min(outerExtNy, nmix)+1);
    if (uy < 1.0) {
      mixingY[jy] = 0.5*(1.0+std::erf(Lmix*(1.0-2.0*uy)/std::sqrt(4.0*uy*(1.0-uy))));
    } else {
      mixingY[jy] = 0.0;
    }
  }

  // Prepare Boyd mask components
  const double Lboyd = params_.Lboyd.value();
  std::vector<double> boydX(outerExtNx);
  for (size_t jx = 0; jx < outerExtNx; ++jx) {
    const double ux = static_cast<double>(jx+1)/static_cast<double>(outerExtNx+1);
    if (ux < 1.0) {
      boydX[jx] = 0.5*(1.0+std::erf(Lboyd*(1.0-2.0*ux)/std::sqrt(4.0*ux*(1.0-ux))));
    } else {
      boydX[jx] = 0.0;
    }
  }
  std::vector<double> boydY(outerExtNy);
  for (size_t jy = 0; jy < outerExtNy; ++jy) {
    const double uy = static_cast<double>(jy+1)/static_cast<double>(outerExtNy+1);
    if (uy < 1.0) {
      boydY[jy] = 0.5*(1.0+std::erf(Lboyd*(1.0-2.0*uy)/std::sqrt(4.0*uy*(1.0-uy))));
    } else {
      boydY[jy] = 0.0;
    }
  }

  // Ghost points
  const auto ghostView = make_view<int, 1>(outerFs.ghost());

  // Index fields views
  const auto indexIView = make_view<int, 1>(outerFs.index_i());
  const auto indexJView = make_view<int, 1>(outerFs.index_j());

  // Loop over local biperiodization operations
  int jx, jy;
  double mixing, boyd;
  localBiperSize_ = 0;
  commBiperSize_ = 0;
  for (int outerJnode = 0; outerJnode < outerFs.size(); ++outerJnode) {
    if (ghostView(outerJnode) == 0) {
      // Grid indices
      jx = indexIView(outerJnode)-1;
      jy = indexJView(outerJnode)-1;

      if (jx < static_cast<int>(physicalNx) && jy < static_cast<int>(physicalNy)) {
        // Physical zone
        addBiperElement(outerJnode, jx, jy, 1.0);
      } else {
        // Indices aliases
        const size_t leftCen = 0;
        const size_t leftSym = outerExtNx-(jx-physicalNx);
        const size_t rightSym = 2*(physicalNx-1)-jx;
        const size_t rightCen = physicalNx-1;
        const size_t bottomCen = 0;
        const size_t bottomSym = outerExtNy-(jy-physicalNy);
        const size_t topCen = physicalNy-1;
        const size_t topSym = 2*(physicalNy-1)-jy;

        // Masks indices
        const size_t left = outerExtNx-(jx-physicalNx)-1;
        const size_t right = jx-physicalNx;
        const size_t bottom = outerExtNy-(jy-physicalNy)-1;
        const size_t top = jy-physicalNy;

        if (jx > static_cast<int>(physicalNx)-1 && jy > static_cast<int>(physicalNy)-1) {
          // Right-top corner

          // Left-bottom corner
          mixing = mixingX[left]*mixingY[bottom];
          boyd = boydX[left]*boydY[bottom];
          addBiperElement(outerJnode, leftCen, bottomCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, leftSym, bottomSym, (1.0-2.0*mixing)*boyd);

          // Left-top corner
          mixing = mixingX[left]*mixingY[top];
          boyd = boydX[left]*boydY[top];
          addBiperElement(outerJnode, leftCen, topCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, leftSym, topSym, (1.0-2.0*mixing)*boyd);

          // Right-top corner
          mixing = mixingX[right]*mixingY[top];
          boyd = boydX[right]*boydY[top];
          addBiperElement(outerJnode, rightCen, topCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, rightSym, topSym, (1.0-2.0*mixing)*boyd);

          // Right-bottom corner
          mixing = mixingX[right]*mixingY[bottom];
          boyd = boydX[right]*boydY[bottom];
          addBiperElement(outerJnode, rightCen, bottomCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, rightSym, bottomSym, (1.0-2.0*mixing)*boyd);
        } else if (jx > static_cast<int>(physicalNx)-1) {
          // Right side

          // Left side
          mixing = mixingX[left];
          boyd = boydX[left];
          addBiperElement(outerJnode, leftCen, jy, 2.0*mixing*boyd);
          addBiperElement(outerJnode, leftSym, jy, (1.0-2.0*mixing)*boyd);

          // Right side
          mixing = mixingX[right];
          boyd = boydX[right];
          addBiperElement(outerJnode, rightCen, jy, 2.0*mixing*boyd);
          addBiperElement(outerJnode, rightSym, jy, (1.0-2.0*mixing)*boyd);
        } else if (jy > static_cast<int>(physicalNy)-1) {
          // Top side

          // Bottom side
          mixing = mixingY[bottom];
          boyd = boydY[bottom];
          addBiperElement(outerJnode, jx, bottomCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, jx, bottomSym, (1.0-2.0*mixing)*boyd);

          // Top side
          mixing = mixingY[top];
          boyd = boydY[top];
          addBiperElement(outerJnode, jx, topCen, 2.0*mixing*boyd);
          addBiperElement(outerJnode, jx, topSym, (1.0-2.0*mixing)*boyd);
        }
      }
    }
  }

  // Reorder according to innerTaskVec_
  std::vector<int> index(commBiperSize_);
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&](int i, int j)
    {return innerTaskVec_[i] < innerTaskVec_[j];});

  // Reorder vectors
  outerJnodeVecOrdered_.resize(commBiperSize_);
  std::vector<size_t> innerJnodeGlbVecOrdered(commBiperSize_);
  std::vector<size_t> innerTaskVecOrdered(commBiperSize_);
  weightVecOrdered_.resize(commBiperSize_);
  for (size_t jj = 0; jj < commBiperSize_; ++jj) {
    outerJnodeVecOrdered_[jj] = outerJnodeVec_[index[jj]];
    innerJnodeGlbVecOrdered[jj] = innerJnodeGlbVec_[index[jj]];
    innerTaskVecOrdered[jj] = innerTaskVec_[index[jj]];
    weightVecOrdered_[jj] = weightVec_[index[jj]];
  }

  // Detect and map duplicated communications
  mappingFull2Red_.resize(commBiperSize_);
  std::vector<int> redInnerJnodeGlb;
  std::vector<int> redInnerTask;
  recvSize_ = 0;
  for (size_t jj = 0; jj < commBiperSize_; ++jj) {
    // Search if element exists
    const auto it = std::find(redInnerJnodeGlb.begin(), redInnerJnodeGlb.end(),
      innerJnodeGlbVecOrdered[jj]);
    if (it == redInnerJnodeGlb.end()) {
      // Add element
      mappingFull2Red_[jj] = recvSize_;
      redInnerJnodeGlb.push_back(innerJnodeGlbVecOrdered[jj]);
      redInnerTask.push_back(innerTaskVecOrdered[jj]);
      ++recvSize_;
    } else {
      // Existing element
      const size_t jjOld = std::distance(redInnerJnodeGlb.begin(), it);
      mappingFull2Red_[jj] = jjOld;
    }
  }

  // RecvCounts
  recvCounts_.resize(comm_.size());
  std::fill(recvCounts_.begin(), recvCounts_.end(), 0);
  for (size_t jjRed = 0; jjRed < recvSize_; ++jjRed) {
    const size_t jt = redInnerTask[jjRed];
    ++recvCounts_[jt];
  }

  // RecvDispls
  recvDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    recvDispls_[jt] = static_cast<int>(jt ? recvDispls_[jt-1] + recvCounts_[jt-1] : 0);
  }

  // SendCounts
  sendCounts_.resize(comm_.size());
  comm_.allToAll(recvCounts_, sendCounts_);

  // SendSize
  sendSize_ = 0;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendSize_ += sendCounts_[jt];
  }

  // SendDispls
  sendDispls_.resize(comm_.size());
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendDispls_[jt] = static_cast<int>(jt ? sendDispls_[jt-1] + sendCounts_[jt-1] : 0);
  }

  // Communicate redInnerJnodeGlb
  std::vector<int> sendInnerJnodeGlb(sendSize_);
  comm_.allToAllv(redInnerJnodeGlb.data(), recvCounts_.data(), recvDispls_.data(),
    sendInnerJnodeGlb.data(), sendCounts_.data(), sendDispls_.data());

  // Get local inner index
  sendInnerJnode_.resize(sendSize_);
  for (size_t jj = 0; jj < sendSize_; ++jj) {
    const int innerGlbJnode = sendInnerJnodeGlb[jj];
    innerGrid_.index2ij(innerGlbJnode, jx,  jy);
    sendInnerJnode_[jj] = innerFs_.index(jx, jy);
  }

  // Scale counts and displs for all levels
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendCounts_[jt] *= nvz_;
    sendDispls_[jt] *= nvz_;
    recvCounts_[jt] *= nvz_;
    recvDispls_[jt] *= nvz_;
  }

  oops::Log::trace() << classname() << "::BiperiodizationImpl done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationImpl::addBiperElement(const size_t & outerJnode,
                                          const size_t & jx,
                                          const size_t & jy,
                                          const double & weight) {
  // Inner global index
  size_t innerJnodeGlb = innerGrid_.index(jx, jy);

  // Inner task
  size_t innerTask = innerPartition_[innerJnodeGlb];

  if (innerTask == myrank_) {
    // Inner local index
    size_t innerJnode = innerFs_.index(jx, jy);

    // Add biperiodization element
    ++localBiperSize_;
    localOuterJnodeVec_.push_back(outerJnode);
    localInnerJnodeVec_.push_back(innerJnode);
    localWeightVec_.push_back(weight);
  } else {
    // Add biperiodization element
    ++commBiperSize_;
    outerJnodeVec_.push_back(outerJnode);
    innerJnodeGlbVec_.push_back(innerJnodeGlb);
    innerTaskVec_.push_back(innerTask);
    weightVec_.push_back(weight);
  }
}

// -----------------------------------------------------------------------------

void BiperiodizationImpl::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Copy inactive fields into output FieldSet
  atlas::FieldSet outputFset;
  for (const auto & field : fset) {
    if (!vars_.has(field.name())) {
      outputFset.add(field);
    }
  }

  // Initialization
  std::vector<double> sendBuf(sendSize_*nvz_);
  size_t zOffset = 0;

  for (const auto & var : vars_) {
    // Check inner field
    const size_t nz = var.getLevels();
    const auto innerField = fset[var.name()];
    ASSERT(innerField.shape(0) == static_cast<int>(innerFs_.size()));
    ASSERT(innerField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto innerView = make_view<double, 2>(innerField);

    // Serialize inner data
    for (size_t jj = 0; jj < sendSize_; ++jj) {
      // Inner local index
      const int innerJnode = sendInnerJnode_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = jj*nvz_ + jvz;

        // Copy data
        sendBuf[jb] = innerView(innerJnode, jz);
      }
    }
    zOffset += nz;
  }

  // Communication
  std::vector<double> recvBuf(recvSize_*nvz_);
  comm_.allToAllv(sendBuf.data(), sendCounts_.data(), sendDispls_.data(),
    recvBuf.data(), recvCounts_.data(), recvDispls_.data());

  // Initialization
  zOffset = 0;

  for (const auto & var : vars_) {
    // Get inner field
    const auto innerField = fset[var.name()];

    // Create outer field
    const size_t nz = var.getLevels();
    atlas::Field outerField = outerGeometryData_.functionSpace().createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    outputFset.add(outerField);

    // Get fields views
    auto outerView = make_view<double, 2>(outerField);
    const auto innerView = make_view<double, 2>(innerField);

    // Deserialize outer data and apply weight
    outerView.assign(0.0);
    for (size_t jj = 0; jj < localBiperSize_; ++jj) {
      // Local indices
      const int innerJnode = localInnerJnodeVec_[jj];
      const int outerJnode = localOuterJnodeVec_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Copy data
        outerView(outerJnode, jz) += innerView(innerJnode, jz)*localWeightVec_[jj];
      }
    }
    for (size_t jj = 0; jj < commBiperSize_; ++jj) {
      // Outer local index
      const int outerJnode = outerJnodeVecOrdered_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = mappingFull2Red_[jj]*nvz_ + jvz;

        // Copy data
        outerView(outerJnode, jz) += recvBuf[jb]*weightVecOrdered_[jj];
      }
    }
    zOffset += nz;
  }

  // Replace FieldSet
  fset = outputFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationImpl::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Copy inactive fields into output FieldSet
  atlas::FieldSet outputFset;
  for (const auto & field : fset) {
    if (!vars_.has(field.name())) {
      outputFset.add(field);
    }
  }

  // Initialization
  std::vector<double> recvBuf(recvSize_*nvz_, 0.0);
  size_t zOffset = 0;

  for (const auto & var : vars_) {
    // Create inner field
    const size_t nz = var.getLevels();
    atlas::Field innerField = innerFs_.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    outputFset.add(innerField);

    // Check outer field
    const auto outerField = fset[var.name()];
    ASSERT(outerField.shape(0) == static_cast<int>(outerGeometryData_.functionSpace().size()));
    ASSERT(outerField.shape(1) == static_cast<int>(nz));

    // Get fields views
    const auto outerView = make_view<double, 2>(outerField);
    auto innerView = make_view<double, 2>(innerField);

    // Deserialize outer data and apply weight
    innerView.assign(0.0);
    for (size_t jj = 0; jj < localBiperSize_; ++jj) {
      // Local indices
      const int innerJnode = localInnerJnodeVec_[jj];
      const int outerJnode = localOuterJnodeVec_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Copy data
        innerView(innerJnode, jz) += outerView(outerJnode, jz)*localWeightVec_[jj];
      }
    }
    for (size_t jj = 0; jj < commBiperSize_; ++jj) {
      // Outer local index
      const int outerJnode = outerJnodeVecOrdered_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = mappingFull2Red_[jj]*nvz_ + jvz;

        // Copy data
        recvBuf[jb] += outerView(outerJnode, jz)*weightVecOrdered_[jj];
      }
    }
    zOffset += nz;
  }

  // Communication
  std::vector<double> sendBuf(sendSize_*nvz_);
  comm_.allToAllv(recvBuf.data(), recvCounts_.data(), recvDispls_.data(),
    sendBuf.data(), sendCounts_.data(), sendDispls_.data());

  // Initialization
  zOffset = 0;

  for (const auto & var : vars_) {
    // Get inner field
    const size_t nz = var.getLevels();
    auto innerField = outputFset[var.name()];

    // Get field view
    auto innerView = make_view<double, 2>(innerField);

    // Serialize inner data
    for (size_t jj = 0; jj < sendSize_; ++jj) {
      // Inner local index
      const int innerJnode = sendInnerJnode_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = jj*nvz_ + jvz;

        // Copy data
        innerView(innerJnode, jz) += sendBuf[jb];
      }
    }
    zOffset += nz;
  }

  // Replace FieldSet
  fset = outputFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationImpl::leftInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  // Copy inactive fields into output FieldSet
  atlas::FieldSet outputFset;
  for (const auto & field : fset) {
    if (!vars_.has(field.name())) {
      outputFset.add(field);
    }
  }

  // Initialization
  std::vector<double> recvBuf(recvSize_*nvz_, 0.0);
  size_t zOffset = 0;

  for (const auto & var : vars_) {
    // Create inner field
    const size_t nz = var.getLevels();
    atlas::Field innerField = innerFs_.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    outputFset.add(innerField);

    // Check outer field
    const auto outerField = fset[var.name()];
    ASSERT(outerField.shape(0) == static_cast<int>(outerGeometryData_.functionSpace().size()));
    ASSERT(outerField.shape(1) == static_cast<int>(nz));

    // Get fields views
    auto innerView = make_view<double, 2>(innerField);
    const auto outerView = make_view<double, 2>(outerField);

    // Deserialize outer data and apply weight
    innerView.assign(0.0);
    for (size_t jj = 0; jj < localBiperSize_; ++jj) {
      // Local indices
      const int innerJnode = localInnerJnodeVec_[jj];
      const int outerJnode = localOuterJnodeVec_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Copy data
        if (oops::is_close_relative(localWeightVec_[jj], 1.0, 1.0e-12, 0,
          oops::TestVerbosity::SILENT)) {
          innerView(innerJnode, jz) += outerView(outerJnode, jz);
        }
      }
    }
    for (size_t jj = 0; jj < commBiperSize_; ++jj) {
      // Outer local index
      const int outerJnode = outerJnodeVecOrdered_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = mappingFull2Red_[jj]*nvz_ + jvz;

        // Copy data
        if (oops::is_close_relative(weightVecOrdered_[jj], 1.0, 1.0e-12, 0,
          oops::TestVerbosity::SILENT)) {
          recvBuf[jb] = outerView(outerJnode, jz);
        }
      }
    }
    zOffset += nz;
  }

  // Communication
  std::vector<double> sendBuf(sendSize_*nvz_);
  comm_.allToAllv(recvBuf.data(), recvCounts_.data(), recvDispls_.data(),
    sendBuf.data(), sendCounts_.data(), sendDispls_.data());

  // Initialization
  zOffset = 0;

  for (const auto & var : vars_) {
    // Get inner field
    const size_t nz = var.getLevels();
    auto innerField = outputFset[var.name()];

    // Get field view
    auto innerView = make_view<double, 2>(innerField);

    // Serialize inner data
    for (size_t jj = 0; jj < sendSize_; ++jj) {
      // Inner local index
      const int innerJnode = sendInnerJnode_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = jj*nvz_ + jvz;

        // Copy data
        innerView(innerJnode, jz) += sendBuf[jb];
      }
    }
    zOffset += nz;
  }

  // Replace FieldSet
  fset = outputFset;

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
