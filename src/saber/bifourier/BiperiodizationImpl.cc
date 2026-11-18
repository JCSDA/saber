/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BiperiodizationImpl.h"

#include <algorithm>

#include "atlas/grid.h"

#include "oops/util/FloatCompare.h"

using atlas::array::make_indexview;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

BiperiodizationData::BiperiodizationData(
       const eckit::mpi::Comm & comm,
       const atlas::functionspace::StructuredColumns & inputFs,
       const std::vector<int> & inputPartition,
       const atlas::functionspace::StructuredColumns & outputFs,
       const size_t & physicalNx,
       const size_t & physicalNy,
       const size_t & outputExtNx,
       const size_t & outputExtNy,
       const size_t & nmix,
       const double & Lmix,
       const double & Lboyd,
       const oops::Variables & vars)
  : comm_(comm),
    inputFs_(inputFs),
    inputPartition_(inputPartition),
    outputFs_(outputFs),
    vars_(vars)
{
  oops::Log::trace() << classname() << "::BiperiodizationData starting" << std::endl;

  // Number of levels for all variables
  nvz_ = 0;
  for (const auto & var : vars_) {
    nvz_ += var.getLevels();
  }

  // Prepare mixing mask components
  const size_t nmixX = nmix == 0 ? outputExtNx : std::min(nmix, outputExtNx);
  std::vector<double> mixingX(outputExtNx);
  for (size_t jx = 0; jx < outputExtNx; ++jx) {
    const double ux = static_cast<double>(jx+1)/static_cast<double>(nmixX+1);
    if (ux < 1.0) {
      mixingX[jx] = 0.5*(1.0+std::erf(Lmix*(1.0-2.0*ux)/std::sqrt(4.0*ux*(1.0-ux))));
    } else {
      mixingX[jx] = 0.0;
    }
  }
  const size_t nmixY = nmix == 0 ? outputExtNy : std::min(nmix, outputExtNy);
  std::vector<double> mixingY(outputExtNy);
  for (size_t jy = 0; jy < outputExtNy; ++jy) {
    const double uy = static_cast<double>(jy+1)/static_cast<double>(std::min(outputExtNy, nmixY)+1);
    if (uy < 1.0) {
      mixingY[jy] = 0.5*(1.0+std::erf(Lmix*(1.0-2.0*uy)/std::sqrt(4.0*uy*(1.0-uy))));
    } else {
      mixingY[jy] = 0.0;
    }
  }

  // Prepare Boyd mask components
  std::vector<double> boydX(outputExtNx);
  for (size_t jx = 0; jx < outputExtNx; ++jx) {
    const double ux = static_cast<double>(jx+1)/static_cast<double>(outputExtNx+1);
    if (ux < 1.0) {
      boydX[jx] = 0.5*(1.0+std::erf(Lboyd*(1.0-2.0*ux)/std::sqrt(4.0*ux*(1.0-ux))));
    } else {
      boydX[jx] = 0.0;
    }
  }
  std::vector<double> boydY(outputExtNy);
  for (size_t jy = 0; jy < outputExtNy; ++jy) {
    const double uy = static_cast<double>(jy+1)/static_cast<double>(outputExtNy+1);
    if (uy < 1.0) {
      boydY[jy] = 0.5*(1.0+std::erf(Lboyd*(1.0-2.0*uy)/std::sqrt(4.0*uy*(1.0-uy))));
    } else {
      boydY[jy] = 0.0;
    }
  }

  // Ghost points
  const auto ghostView = make_view<int, 1>(outputFs_.ghost());

  // Index fields views
  const auto indexIView = make_indexview<int, 1>(outputFs_.index_i());
  const auto indexJView = make_indexview<int, 1>(outputFs_.index_j());

  // Loop over local biperiodization operations
  int jx, jy;
  double mixing, boyd;
  localBiperSize_ = 0;
  commBiperSize_ = 0;
  for (int outputJnode = 0; outputJnode < outputFs_.size(); ++outputJnode) {
    if (ghostView(outputJnode) == 0) {
      // Grid indices
      jx = indexIView(outputJnode);
      jy = indexJView(outputJnode);

      // Copy outputJnode
      outputJnode_ = outputJnode;

      if (jx < static_cast<int>(physicalNx) && jy < static_cast<int>(physicalNy)) {
        // Physical zone
        addBiperElement(jx, jy, 1.0);
      } else {
        // Indices aliases
        const size_t leftCen = 0;
        const size_t leftSym = outputExtNx-(jx-physicalNx);
        const size_t rightSym = 2*(physicalNx-1)-jx;
        const size_t rightCen = physicalNx-1;
        const size_t bottomCen = 0;
        const size_t bottomSym = outputExtNy-(jy-physicalNy);
        const size_t topCen = physicalNy-1;
        const size_t topSym = 2*(physicalNy-1)-jy;

        // Masks indices
        const size_t left = outputExtNx-(jx-physicalNx)-1;
        const size_t right = jx-physicalNx;
        const size_t bottom = outputExtNy-(jy-physicalNy)-1;
        const size_t top = jy-physicalNy;

        if (jx > static_cast<int>(physicalNx)-1 && jy > static_cast<int>(physicalNy)-1) {
          // Right-top corner

          // Left-bottom corner
          mixing = mixingX[left]*mixingY[bottom];
          boyd = boydX[left]*boydY[bottom];
          addBiperElement(leftCen, bottomCen, 2.0*mixing*boyd);
          addBiperElement(leftSym, bottomSym, (1.0-2.0*mixing)*boyd);

          // Left-top corner
          mixing = mixingX[left]*mixingY[top];
          boyd = boydX[left]*boydY[top];
          addBiperElement(leftCen, topCen, 2.0*mixing*boyd);
          addBiperElement(leftSym, topSym, (1.0-2.0*mixing)*boyd);

          // Right-top corner
          mixing = mixingX[right]*mixingY[top];
          boyd = boydX[right]*boydY[top];
          addBiperElement(rightCen, topCen, 2.0*mixing*boyd);
          addBiperElement(rightSym, topSym, (1.0-2.0*mixing)*boyd);

          // Right-bottom corner
          mixing = mixingX[right]*mixingY[bottom];
          boyd = boydX[right]*boydY[bottom];
          addBiperElement(rightCen, bottomCen, 2.0*mixing*boyd);
          addBiperElement(rightSym, bottomSym, (1.0-2.0*mixing)*boyd);
        } else if (jx > static_cast<int>(physicalNx)-1) {
          // Right side

          // Left side
          mixing = mixingX[left];
          boyd = boydX[left];
          addBiperElement(leftCen, jy, 2.0*mixing*boyd);
          addBiperElement(leftSym, jy, (1.0-2.0*mixing)*boyd);

          // Right side
          mixing = mixingX[right];
          boyd = boydX[right];
          addBiperElement(rightCen, jy, 2.0*mixing*boyd);
          addBiperElement(rightSym, jy, (1.0-2.0*mixing)*boyd);
        } else if (jy > static_cast<int>(physicalNy)-1) {
          // Top side

          // Bottom side
          mixing = mixingY[bottom];
          boyd = boydY[bottom];
          addBiperElement(jx, bottomCen, 2.0*mixing*boyd);
          addBiperElement(jx, bottomSym, (1.0-2.0*mixing)*boyd);

          // Top side
          mixing = mixingY[top];
          boyd = boydY[top];
          addBiperElement(jx, topCen, 2.0*mixing*boyd);
          addBiperElement(jx, topSym, (1.0-2.0*mixing)*boyd);
        }
      }
    }
  }

  // Reorder according to inputTaskVec
  std::vector<int> index(commBiperSize_);
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&](int i, int j)
    {return inputTaskVec_[i] < inputTaskVec_[j];});

  // Reorder vectors
  outputJnodeVecOrdered_.resize(commBiperSize_);
  std::vector<size_t> inputJnodeGlbVecOrdered(commBiperSize_);
  std::vector<size_t> inputTaskVecOrdered(commBiperSize_);
  weightVecOrdered_.resize(commBiperSize_);
  for (size_t jj = 0; jj < commBiperSize_; ++jj) {
    outputJnodeVecOrdered_[jj] = outputJnodeVec_[index[jj]];
    inputJnodeGlbVecOrdered[jj] = inputJnodeGlbVec_[index[jj]];
    inputTaskVecOrdered[jj] = inputTaskVec_[index[jj]];
    weightVecOrdered_[jj] = weightVec_[index[jj]];
  }

  // Detect and map duplicated communications
  mappingFull2Red_.resize(commBiperSize_);
  std::vector<int> redInputJnodeGlb;
  std::vector<int> redInputTask;
  recvSize_ = 0;
  for (size_t jj = 0; jj < commBiperSize_; ++jj) {
    // Search if element exists
    const auto it = std::find(redInputJnodeGlb.begin(), redInputJnodeGlb.end(),
      inputJnodeGlbVecOrdered[jj]);
    if (it == redInputJnodeGlb.end()) {
      // Add element
      mappingFull2Red_[jj] = recvSize_;
      redInputJnodeGlb.push_back(inputJnodeGlbVecOrdered[jj]);
      redInputTask.push_back(inputTaskVecOrdered[jj]);
      ++recvSize_;
    } else {
      // Existing element
      const size_t jjOld = std::distance(redInputJnodeGlb.begin(), it);
      mappingFull2Red_[jj] = jjOld;
    }
  }

  // RecvCounts
  recvCounts_.resize(comm_.size(), 0);
  for (size_t jjRed = 0; jjRed < recvSize_; ++jjRed) {
    const size_t jt = redInputTask[jjRed];
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

  // Communicate redInputJnodeGlb
  std::vector<int> sendInputJnodeGlb(sendSize_);
  comm_.allToAllv(redInputJnodeGlb.data(), recvCounts_.data(), recvDispls_.data(),
    sendInputJnodeGlb.data(), sendCounts_.data(), sendDispls_.data());

  // Get local input index
  sendInputJnode_.resize(sendSize_);
  for (size_t jj = 0; jj < sendSize_; ++jj) {
    const int inputGlbJnode = sendInputJnodeGlb[jj];
    inputFs_.grid().index2ij(inputGlbJnode, jx,  jy);
    sendInputJnode_[jj] = inputFs_.index(jx, jy);
  }

  // Scale counts and displs for all levels
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    sendCounts_[jt] *= nvz_;
    sendDispls_[jt] *= nvz_;
    recvCounts_[jt] *= nvz_;
    recvDispls_[jt] *= nvz_;
  }

  oops::Log::trace() << classname() << "::BiperiodizationData done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationData::multiply(atlas::FieldSet & fset) const {
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
    // Check input field
    const size_t nz = var.getLevels();
    const auto inputField = fset[var.name()];
    ASSERT(inputField.shape(0) == static_cast<int>(inputFs_.size()));
    ASSERT(inputField.shape(1) == static_cast<int>(nz));

    // Get field view
    const auto inputView = make_view<double, 2>(inputField);

    // Serialize input data
    for (size_t jj = 0; jj < sendSize_; ++jj) {
      // Input local index
      const int inputJnode = sendInputJnode_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = jj*nvz_ + jvz;

        // Copy data
        sendBuf[jb] = inputView(inputJnode, jz);
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
    // Get input field
    const auto inputField = fset[var.name()];

    // Create output field
    const size_t nz = var.getLevels();
    atlas::Field outputField = outputFs_.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    outputFset.add(outputField);

    // Get fields views
    auto outputView = make_view<double, 2>(outputField);
    const auto inputView = make_view<double, 2>(inputField);

    // Deserialize output data and apply weight
    outputView.assign(0.0);
    for (size_t jj = 0; jj < localBiperSize_; ++jj) {
      // Local indices
      const int inputJnode = localInputJnodeVec_[jj];
      const int outputJnode = localOutputJnodeVec_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Copy data
        outputView(outputJnode, jz) += inputView(inputJnode, jz)*localWeightVec_[jj];
      }
    }
    for (size_t jj = 0; jj < commBiperSize_; ++jj) {
      // Output local index
      const int outputJnode = outputJnodeVecOrdered_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = mappingFull2Red_[jj]*nvz_ + jvz;

        // Copy data
        outputView(outputJnode, jz) += recvBuf[jb]*weightVecOrdered_[jj];
      }
    }
    zOffset += nz;
  }

  // Replace FieldSet
  fset = outputFset;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationData::multiplyAD(atlas::FieldSet & fset) const {
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
    // Create input field
    const size_t nz = var.getLevels();
    atlas::Field inputField = inputFs_.createField<double>(
      atlas::option::name(var.name()) | atlas::option::levels(nz));
    outputFset.add(inputField);

    // Check output field
    const auto outputField = fset[var.name()];
    ASSERT(outputField.shape(0) == static_cast<int>(outputFs_.size()));
    ASSERT(outputField.shape(1) == static_cast<int>(nz));

    // Get fields views
    const auto outputView = make_view<double, 2>(outputField);
    auto inputView = make_view<double, 2>(inputField);

    // Deserialize output data and apply weight
    inputView.assign(0.0);
    for (size_t jj = 0; jj < localBiperSize_; ++jj) {
      // Local indices
      const int inputJnode = localInputJnodeVec_[jj];
      const int outputJnode = localOutputJnodeVec_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Copy data
        inputView(inputJnode, jz) += outputView(outputJnode, jz)*localWeightVec_[jj];
      }
    }
    for (size_t jj = 0; jj < commBiperSize_; ++jj) {
      // Output local index
      const int outputJnode = outputJnodeVecOrdered_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = mappingFull2Red_[jj]*nvz_ + jvz;

        // Copy data
        recvBuf[jb] += outputView(outputJnode, jz)*weightVecOrdered_[jj];
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
    // Get input field
    const size_t nz = var.getLevels();
    auto inputField = outputFset[var.name()];

    // Get field view
    auto inputView = make_view<double, 2>(inputField);

    // Serialize input data
    for (size_t jj = 0; jj < sendSize_; ++jj) {
      // Input local index
      const int inputJnode = sendInputJnode_[jj];

      for (size_t jz = 0; jz < nz; ++jz) {
        // Total level index
        const size_t jvz = zOffset + jz;

        // Buffer index
        size_t jb = jj*nvz_ + jvz;

        // Copy data
        inputView(inputJnode, jz) += sendBuf[jb];
      }
    }
    zOffset += nz;
  }

  // Replace FieldSet
  fset = outputFset;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BiperiodizationData::addBiperElement(const size_t & jx,
                                          const size_t & jy,
                                          const double & weight) {
  // Input global index
  size_t inputJnodeGlb = inputFs_.grid().index(jx, jy);

  // Input task
  size_t inputTask = inputPartition_[inputJnodeGlb];

  if (inputTask == comm_.rank()) {
    // Input local index
    size_t inputJnode = inputFs_.index(jx, jy);

    // Add biperiodization element
    ++localBiperSize_;
    localOutputJnodeVec_.push_back(outputJnode_);
    localInputJnodeVec_.push_back(inputJnode);
    localWeightVec_.push_back(weight);
  } else {
    // Add biperiodization element
    ++commBiperSize_;
    outputJnodeVec_.push_back(outputJnode_);
    inputJnodeGlbVec_.push_back(inputJnodeGlb);
    inputTaskVec_.push_back(inputTask);
    weightVec_.push_back(weight);
  }
}

// -----------------------------------------------------------------------------

BiperiodizationImpl::BiperiodizationImpl(const oops::GeometryData & outerGeometryData,
                                         const oops::Variables & vars,
                                         const Parameters_ & params)
  : comm_(outerGeometryData.comm())
{
  oops::Log::trace() << classname() << "::BiperiodizationImpl starting" << std::endl;

  // Check function space type
  ASSERT(outerGeometryData.functionSpace().type() == "StructuredColumns");

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
  const size_t innerExtNx = params.innerExtNx.value();
  const size_t innerExtNy = params.innerExtNy.value();
  const size_t outerExtNx = params.outerExtNx.value();
  const size_t outerExtNy = params.outerExtNy.value();

  // Define physical grid size
  const size_t physicalNx = outerNx - outerExtNx;
  const size_t physicalNy = outerNy - outerExtNy;

  // Allocate outer partition
  std::vector<int> outerPartition(outerGrid.size());

  // Communicate partition
  atlas::Field glbPartitionField = outerFs.createField<int>(
    atlas::option::name("global partition") | atlas::option::global());
  outerFs.gather(outerFs.partition(), glbPartitionField);

  if (comm_.rank() == 0) {
    // Fill partition vector
    const auto glbPartitionView = make_view<int, 1>(glbPartitionField);
    for (int jnodeGlb = 0; jnodeGlb < outerGrid.size(); ++jnodeGlb) {
      outerPartition[jnodeGlb] = glbPartitionView(jnodeGlb);
    }
  }

  // Broadcast partition vector
  comm_.broadcast(outerPartition.begin(), outerPartition.end(), 0);

  // Create inner grid
  atlas::StructuredGrid innerGrid;

  // Create inner partition
  std::vector<int> innerPartition;

  if (innerExtNx == outerExtNx && innerExtNy == outerExtNy) {
    // Same function space
    oops::Log::info() << "Info     : Inner grid resolution = outer grid resolution" << std::endl;
    sameFs_ = true;

    // Copy grid
    innerGrid = outerGrid;

    // Copy function space
    innerFs_ = outerFs;

    // Allocate inner partition
    innerPartition.resize(innerGrid.size());

    // Copy outer partition
    innerPartition = outerPartition;
  } else {
    // Same function space
    sameFs_ = false;

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
    innerGrid = atlas::StructuredGrid(innerXSpace, innerYSpace, outerFs.grid().projection());
    oops::Log::test() << "Info     : Inner grid size: " << innerGrid.size() << std::endl;

    // Get outer Partitioner name
    const std::string outerPartitionerName = outerFs.distribution();
    oops::Log::info() << "Info     : Outer partitioner: " << outerPartitionerName << std::endl;

    // Get inner Partitioner name
    std::string innerPartitionerName;
    if (outerPartitionerName == "custom") {
      // Mandatory input parameter (custom is not a valid inner partitioner name)
      ASSERT(params.innerPartitioner.value());
      innerPartitionerName = *params.innerPartitioner.value();
    } else {
      // Optional input parameter
      if (params.innerPartitioner.value()) {
        innerPartitionerName = *params.innerPartitioner.value();
      } else {
        innerPartitionerName = outerPartitionerName;
      }
    }

    // Get inner Partitioner
    const auto innerPartitioner =  atlas::grid::Partitioner(innerPartitionerName);

    // Allocate inner partition
    innerPartition.resize(innerGrid.size());

    // Generate inner partition
    innerPartitioner.partition(innerGrid, innerPartition.data());

    // Generate inner Distribution
    atlas::grid::Distribution innerDistribution(comm_.size(), innerGrid.size(),
    innerPartition.data());

    // Generate inner FunctionSpace
    innerFs_ = atlas::functionspace::StructuredColumns(innerGrid, innerDistribution);
  }

  // Prepare mask components
  const size_t nmix = params.nmix.value();
  const double Lmix = params.Lmix.value();
  const double Lboyd = params.Lboyd.value();

  // Prepare inner to outer
  direct_ = std::make_unique<BiperiodizationData>(comm_, innerFs_, innerPartition, outerFs,
    physicalNx, physicalNy, outerExtNx, outerExtNy, nmix, Lmix, Lboyd, vars);

  // Prepare outer to inner
  inverse_ = std::make_unique<BiperiodizationData>(comm_, outerFs, outerPartition, innerFs_,
    physicalNx, physicalNy, innerExtNx, innerExtNy, nmix, Lmix, Lboyd, vars);

  oops::Log::trace() << classname() << "::BiperiodizationImpl done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
