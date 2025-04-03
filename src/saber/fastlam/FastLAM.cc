/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/FastLAM.h"

#include <math.h>
#include <netcdf.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/RandomField.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<FastLAM> makerFastLAM_("FastLAM");

// -----------------------------------------------------------------------------

FastLAM::FastLAM(const oops::GeometryData & gdata,
                 const oops::Variables & activeVars,
                 const eckit::Configuration & covarConf,
                 const Parameters_ & params,
                 const oops::FieldSet3D & xb,
                 const oops::FieldSet3D & fg) :
    SaberCentralBlockBase(params, xb.validTime()),
    gdata_(gdata),
    comm_(gdata_.comm()),
    activeVars_(activeVars),
    params_(params.calibration.value() != boost::none ? *params.calibration.value()
      : *params.read.value()),
    fieldsMetaData_(params.fieldsMetaData.value())
{
  oops::Log::trace() << classname() << "::FastLAM starting" << std::endl;
  // Check function space type
  ASSERT(gdata_.functionSpace().type() == "StructuredColumns");

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Index fields
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  const auto indexX0View = atlas::array::make_view<int, 1>(fs.index_i());
  const auto indexY0View = atlas::array::make_view<int, 1>(fs.index_j());

  // Get grid size
  nx0_ = 0;
  ny0_ = 0;
  nodes0_ = fs.size();
  for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
    if (ghostView(jnode0) == 0) {
      nx0_ = std::max(nx0_, static_cast<size_t>(indexX0View(jnode0)));
      ny0_ = std::max(ny0_, static_cast<size_t>(indexY0View(jnode0)));
    }
  }
  comm_.allReduceInPlace(nx0_, eckit::mpi::max());
  comm_.allReduceInPlace(ny0_, eckit::mpi::max());
  oops::Log::info() << "Info     : Regional grid size: " << nx0_ << "x" << ny0_ << std::endl;

  // Define 2d active variables
  active2dVars_ = oops::Variables();
  for (const auto & var : activeVars_) {
    if (var.getLevels() == 1) {
      active2dVars_.push_back(var);
    }
  }

  // Create groups
  if (params_.groups.value() == boost::none) {
    // No group specified, each variable is its own group
    for (const auto & var : activeVars_) {
      // Define group properties
      Group group;
      // TODO(AS): I think name_, nz0_ can be removed if variables_ are oops::Variables
      group.name_ = var.name();
      group.nz0_ = var.getLevels();
      group.varInModelFile_ = var.name();
      group.variables_ = {var.name()};

      // Add group
      groups_.push_back(group);
    }
  } else {
    // Copy groups
    for (const auto & groupParams : *params_.groups.value()) {
      // Define group properties
      Group group;
      group.name_ = groupParams.name.value();
      group.nz0_ = 1;
      for (const auto & var : groupParams.variables.value()) {
        if (!active2dVars_.has(var)) {
          if (group.nz0_ == 1) {
            // Assign number of levels
            group.nz0_ = static_cast<size_t>(activeVars_[var].getLevels());
          } else {
            // Check number of levels
            ASSERT(static_cast<int>(group.nz0_) == activeVars_[var].getLevels());
          }
        }
      }
      group.varInModelFile_ = groupParams.varInModelFile.value().get_value_or(group.name_);
      group.variables_ = groupParams.variables.value();

      // Add group
      groups_.push_back(group);
    }
  }

  oops::Log::trace() << classname() << "::FastLAM done" << std::endl;
}

// -----------------------------------------------------------------------------

FastLAM::~FastLAM() {
  oops::Log::trace() << classname() << "::~FastLAM starting" << std::endl;
  oops::Log::trace() << classname() << "::~FastLAM done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Create control vector
  atlas::Field cv("genericCtlVec", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(ctlVecSize()));

  // Sizes, sendcounts and displs
  std::vector<int> sendcounts(comm_.size());
  comm_.allGather(static_cast<int>(ctlVecSize()), sendcounts.begin(), sendcounts.end());
  size_t ctlVecSizeGlb = 0;
  for (const auto ctlVecSize : sendcounts) {
    ctlVecSizeGlb += ctlVecSize;
  }
  std::vector<int> displs(comm_.size());
  displs[0] = 0;
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    displs[jt+1] = displs[jt]+sendcounts[jt];
  }

  // Generate global random vector
  std::vector<double> rand_vec_glb;
  if (comm_.rank() == 0) {
    util::NormalDistributionField dist(ctlVecSizeGlb, 0.0, 1.0);
    rand_vec_glb.resize(ctlVecSizeGlb);
    for (size_t jcv = 0; jcv < ctlVecSizeGlb; ++jcv) {
      rand_vec_glb[jcv] = dist[jcv];
    }
  }

  // Scatter random vector
  std::vector<double> rand_vec(ctlVecSize());
  comm_.scatterv(rand_vec_glb.begin(), rand_vec_glb.end(), sendcounts, displs,
    rand_vec.begin(), rand_vec.end(), 0);

  // Fill control vector
  auto cvView = atlas::array::make_view<double, 1>(cv);
  for (size_t jcv = 0; jcv < ctlVecSize(); ++jcv) {
    cvView(jcv) = rand_vec[jcv];
  }

  // Square-root multiplication
  const size_t index = 0;
  multiplySqrt(cv, fset, index);

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Create control vector
  atlas::Field cv("genericCtlVec", atlas::array::make_datatype<double>(),
    atlas::array::make_shape(ctlVecSize()));
  const size_t index = 0;

  // Square-root multiplication, adjoint
  multiplySqrtAD(fset, cv, index);

  // Square-root multiplication
  multiplySqrt(cv, fset, index);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t FastLAM::ctlVecSize() const {
  oops::Log::trace() << classname() << "::ctlVecSize starting" << std::endl;

  // Loop over bins
  size_t ctlVecSize = 0;
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    // Loop over groups
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      if (params_.strategy.value() == "univariate") {
        // Univariate strategy
        ctlVecSize += data_[jg][jBin]->ctlVecSize()*groups_[jg].variables_.size();
      } else if (params_.strategy.value() == "duplicated") {
        // Duplicated strategy
        ctlVecSize += data_[jg][jBin]->ctlVecSize();
      } else if (params_.strategy.value() == "crossed") {
        // Crossed strategy
        if (jg == 0) {
          ctlVecSize += data_[jg][jBin]->ctlVecSize();
        } else {
          ASSERT(data_[jg][jBin]->ctlVecSize() == data_[0][jBin]->ctlVecSize());
        }
      } else {
        // Wrong multivariate strategy
        throw eckit::UserError("wrong multivariate strategy: " + params_.strategy.value(), Here());
      }
    }
  }

  oops::Log::trace() << classname() << "::ctlVecSize done" << std::endl;
  return ctlVecSize;
}

// -----------------------------------------------------------------------------

void FastLAM::multiplySqrt(const atlas::Field & cv,
                           oops::FieldSet3D & fset,
                           const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Save input FieldSet
  oops::FieldSet3D fsetIn(fset);
  fset.zero();

  // Initialize control vector index
  int index = offset;

  // Loop over bins
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    // Copy input FieldSet
    oops::FieldSet3D fsetBin(fsetIn);

    if (params_.strategy.value() == "univariate") {
      // Univariate strategy
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Loop over variables
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Layer multiplication
          atlas::Field binField = fsetBin[var];
          data_[jg][jBin]->multiplySqrt(cv, binField, index);

          // Update control vector index
          index += data_[jg][jBin]->ctlVecSize();

          // Apply weight square-root and normalization
          auto binView = atlas::array::make_view<double, 2>(binField);
          const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
          const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
          const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
          const auto normView = atlas::array::make_view<double, 2>(normField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                binView(jnode0, jz0) *= wgtSqrtView(jnode0, z0Offset+jz0)
                  *normView(jnode0, z0Offset+jz0);
              }
            }
          }
        }
      }
    } else if (params_.strategy.value() == "duplicated") {
      // Duplicated strategy
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Create group field
        atlas::Field grpField = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto grpView = atlas::array::make_view<double, 2>(grpField);

        // Layer multiplication
        data_[jg][jBin]->multiplySqrt(cv, grpField, index);

        // Update control vector index
        index += data_[jg][jBin]->ctlVecSize();

        // Apply weight square-root and normalization
        const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
        const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
        const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
        const auto normView = atlas::array::make_view<double, 2>(normField);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (ghostView(jnode0) == 0) {
            for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
              grpView(jnode0, jz0) *= wgtSqrtView(jnode0, jz0)*normView(jnode0, jz0);
            }
          }
        }

        // Copy result on all variables of the group
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Copy group field
          atlas::Field binField = fsetBin[var];
          auto binView = atlas::array::make_view<double, 2>(binField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                binView(jnode0, jz0) = grpView(jnode0, z0Offset+jz0);
              }
            }
          }
        }
      }
    } else if (params_.strategy.value() == "crossed") {
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Create group field
        atlas::Field grpField = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto grpView = atlas::array::make_view<double, 2>(grpField);

        // Layer square-root multiplication
        data_[jg][jBin]->multiplySqrt(cv, grpField, index);

        // Apply weight square-root and normalization
        const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
        const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
        const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
        const auto normView = atlas::array::make_view<double, 2>(normField);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (ghostView(jnode0) == 0) {
            for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
              grpView(jnode0, jz0) *= wgtSqrtView(jnode0, jz0)*normView(jnode0, jz0);
            }
          }
        }

        // Copy result on all variables of the group
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Copy group field
          atlas::Field binField = fsetBin[var];
          auto binView = atlas::array::make_view<double, 2>(binField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                binView(jnode0, jz0) = grpView(jnode0, z0Offset+jz0);
              }
            }
          }
        }
      }

      // Update control vector index
      index += data_[0][jBin]->ctlVecSize();
    } else {
      // Wrong multivariate strategy
      throw eckit::UserError("wrong multivariate strategy: " + params_.strategy.value(), Here());
    }

    // Add component
    fset += fsetBin;
  }

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::multiplySqrtAD(const oops::FieldSet3D & fset,
                             atlas::Field & cv,
                             const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Save input FieldSet
  oops::FieldSet3D fsetIn(fset);

  if (params_.strategy.value() == "crossed") {
    // Initialize control vector
    auto cvView = atlas::array::make_view<double, 1>(cv);
    cvView.assign(0.0);
  }

  // Initialize control vector index
  int index = offset;

  // Loop over bins
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    // Copy input FieldSet
    oops::FieldSet3D fsetBin(fsetIn);

    if (params_.strategy.value() == "univariate") {
      // Univariate strategy
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Loop over variables
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Apply weight square-root and normalization
          atlas::Field binField = fsetBin[var];
          auto binView = atlas::array::make_view<double, 2>(binField);
          const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
          const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
          const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
          const auto normView = atlas::array::make_view<double, 2>(normField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                binView(jnode0, jz0) *= wgtSqrtView(jnode0, z0Offset+jz0)
                  *normView(jnode0, z0Offset+jz0);
              }
            }
          }

          // Layer multiplication
          data_[jg][jBin]->multiplySqrtTrans(binField, cv, index);

          // Update control vector index
          index += data_[jg][jBin]->ctlVecSize();
        }
      }
    } else if (params_.strategy.value() == "duplicated") {
      // Duplicated strategy
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Create group field
        atlas::Field grpField = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto grpView = atlas::array::make_view<double, 2>(grpField);
        grpView.assign(0.0);

        // Sum all variables of the group
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Add variable field
          const atlas::Field binField = fsetBin[var];
          const auto binView = atlas::array::make_view<double, 2>(binField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                grpView(jnode0, z0Offset+jz0) += binView(jnode0, jz0);
              }
            }
          }
        }

        // Apply weight square-root and normalization
        const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
        const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
        const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
        const auto normView = atlas::array::make_view<double, 2>(normField);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (ghostView(jnode0) == 0) {
            for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
              grpView(jnode0, jz0) *= wgtSqrtView(jnode0, jz0)*normView(jnode0, jz0);
            }
          }
        }

        // Layer multiplication
        data_[jg][jBin]->multiplySqrtTrans(grpField, cv, index);

        // Update control vector index
        index += data_[jg][jBin]->ctlVecSize();
      }
    } else if (params_.strategy.value() == "crossed") {
      // Crossed strategy
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Initialize temporary control vector
        atlas::Field cvBin("genericCtlVecBin", atlas::array::make_datatype<double>(),
          atlas::array::make_shape(data_[jg][jBin]->ctlVecSize()));
        auto cvBinView = atlas::array::make_view<double, 1>(cvBin);
        ASSERT(data_[jg][jBin]->ctlVecSize() == data_[0][jBin]->ctlVecSize());

        // Create group field
        atlas::Field grpField = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto grpView = atlas::array::make_view<double, 2>(grpField);
        grpView.assign(0.0);

        // Sum all variables of the group
        for (const auto & var : groups_[jg].variables_) {
          // Variable properties
          const size_t varNz0 = activeVars_[var].getLevels();
          const size_t z0Offset = getZ0Offset(var);

          // Add variable field
          const atlas::Field binField = fsetBin[var];
          const auto binView = atlas::array::make_view<double, 2>(binField);
          for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
            if (ghostView(jnode0) == 0) {
              for (size_t jz0 = 0; jz0 < varNz0; ++jz0) {
                grpView(jnode0, z0Offset+jz0) += binView(jnode0, jz0);
              }
            }
          }
        }

        // Apply weight square-root and normalization
        const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
        const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
        const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
        const auto normView = atlas::array::make_view<double, 2>(normField);
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (ghostView(jnode0) == 0) {
            for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
              grpView(jnode0, jz0) *= wgtSqrtView(jnode0, jz0)*normView(jnode0, jz0);
            }
          }
        }

        // Layer square-root multiplication, adjoint
        data_[jg][jBin]->multiplySqrtTrans(grpField, cvBin, 0);

        // Add contribution
        auto cvView = atlas::array::make_view<double, 1>(cv);
        for (size_t jcv = 0; jcv < data_[jg][jBin]->ctlVecSize(); ++jcv) {
          cvView(index+jcv) += cvBinView(jcv);
        }
      }

      // Update control vector index
      index += data_[0][jBin]->ctlVecSize();
    } else {
      // Wrong multivariate strategy
      throw eckit::UserError("wrong multivariate strategy: " + params_.strategy.value(), Here());
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> FastLAM::getReadConfs() const {
  oops::Log::trace() << classname() << "::getReadConfs starting" << std::endl;

  std::vector<std::pair<std::string, eckit::LocalConfiguration>> inputs;
  if (params_.inputModelFilesConf.value() != boost::none) {
    for (const auto & conf : *params_.inputModelFilesConf.value()) {
      // Get parameter
      const std::string param = conf.getString("parameter");

      // Number of components
      size_t nWeight = 0;
      size_t nNormalization = 0;

      if (param == "weight") {
        // Specific case for weight
        nWeight = conf.getInt("number of components");
        for (size_t jBin = 0; jBin < nWeight; ++jBin) {
          // File name
          const std::string fileName = param + "_" + std::to_string(jBin);

          // Get file configuration
          eckit::LocalConfiguration fileConf = getFileConf(comm_, conf);

          // Update configuration
          util::seekAndReplace(fileConf, "%component%", jBin, 2);

          // Add pair
          inputs.push_back(std::make_pair(fileName, fileConf));
        }
      } else if (param == "normalization") {
        // Specific case for normalization
        nNormalization = conf.getInt("number of components");
        for (size_t jBin = 0; jBin < nNormalization; ++jBin) {
          // File name
          const std::string fileName = param + "_" + std::to_string(jBin);

          // Get file configuration
          eckit::LocalConfiguration fileConf = getFileConf(comm_, conf);

          // Update configuration
          util::seekAndReplace(fileConf, "%component%", jBin, 2);

          // Add pair
          inputs.push_back(std::make_pair(fileName, fileConf));
        }
      } else {
        // File name
        const std::string fileName = param;

        // Get file configuration
        eckit::LocalConfiguration fileConf = getFileConf(comm_, conf);

        // Add pair
        inputs.push_back(std::make_pair(fileName, fileConf));
      }
    }
  }

  oops::Log::trace() << classname() << "::getReadConfs done" << std::endl;
  return inputs;
}

// -----------------------------------------------------------------------------

void FastLAM::setReadFields(const std::vector<oops::FieldSet3D> & fsetVec) {
  oops::Log::trace() << classname() << "::setReadFields starting" << std::endl;

  // Get number of weight and normalization components from input files
  size_t nWeight = 0;
  size_t nLayers = 0;
  for (size_t ji = 0; ji < fsetVec.size(); ++ji) {
    for (size_t jBin = 0; jBin < 100; ++jBin) {
      const std::string weightName = "weight_" + std::to_string(jBin);
      if (fsetVec[ji].name() == weightName) {
        ++nWeight;
      }
      const std::string normalizationName = "normalization_" + std::to_string(jBin);
      if (fsetVec[ji].name() == normalizationName) {
        ++nLayers;
      }
    }
  }

  if (nLayers > 0) {
    // Check consistency
    if (nWeight == 0) {
      ASSERT(nLayers == 1);
    } else {
      ASSERT(nLayers == nWeight);
    }

    // Allocate data
    data_.resize(groups_.size());
    weight_.resize(nLayers);
    normalization_.resize(nLayers);
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      // Create layers
      for (size_t jBin = 0; jBin < nLayers; ++jBin) {
        data_[jg].emplace_back(LayerFactory::create(params_, fieldsMetaData_, gdata_,
          groups_[jg].name_, groups_[jg].variables_, nx0_, ny0_, groups_[jg].nz0_));
      }
    }

    for (size_t ji = 0; ji < fsetVec.size(); ++ji) {
      for (size_t jBin = 0; jBin < 100; ++jBin) {
        // Get weight from input files
        std::string weightName = "weight_" + std::to_string(jBin);
        if (fsetVec[ji].name() == weightName) {
          weight_[jBin].reset(new oops::FieldSet3D(validTime_, comm_));
          for (size_t jg = 0; jg < groups_.size(); ++jg) {
            // Copy field
            atlas::Field field = gdata_.functionSpace().createField<double>(
              atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
            auto view = atlas::array::make_view<double, 2>(field);
            atlas::Field inputField = fsetVec[ji][groups_[jg].varInModelFile_];
            auto inputView = atlas::array::make_view<double, 2>(inputField);
            view.assign(inputView);
            weight_[jBin]->add(field);
          }
        }

        // Get normalization components from input files
        std::string normalizationName = "normalization_" + std::to_string(jBin);
        if (fsetVec[ji].name() == normalizationName) {
          normalization_[jBin].reset(new oops::FieldSet3D(validTime_, comm_));
          for (size_t jg = 0; jg < groups_.size(); ++jg) {
            // Copy field
            atlas::Field field = gdata_.functionSpace().createField<double>(
              atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
            auto view = atlas::array::make_view<double, 2>(field);
            atlas::Field inputField = fsetVec[ji][groups_[jg].varInModelFile_];
            auto inputView = atlas::array::make_view<double, 2>(inputField);
            view.assign(inputView);
            normalization_[jBin]->add(field);
          }
        }
      }
    }

    // Set weight to one if not present
    if (nWeight == 0) {
      weight_[0].reset(new oops::FieldSet3D(validTime_, comm_));
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Copy field
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto view = atlas::array::make_view<double, 2>(field);
        view.assign(1.0);
        weight_[0]->add(field);
      }
    }
  }

  for (size_t ji = 0; ji < fsetVec.size(); ++ji) {
    // Get rh from input files
    if (fsetVec[ji].name() == "rh") {
      rh_.reset(new oops::FieldSet3D(validTime_, comm_));
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Copy field
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto view = atlas::array::make_view<double, 2>(field);
        if (fsetVec[ji].has(groups_[jg].varInModelFile_)) {
          const atlas::Field inputField = fsetVec[ji][groups_[jg].varInModelFile_];
          const auto inputView = atlas::array::make_view<double, 2>(inputField);
          view.assign(inputView);
        } else {
          throw eckit::Exception("missing " + groups_[jg].varInModelFile_ + " in model file",
            Here());
        }
        rh_->add(field);
      }
    }

    // Get rv from input files
    if (fsetVec[ji].name() == "rv") {
      rv_.reset(new oops::FieldSet3D(validTime_, comm_));
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        // Copy field
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
        auto view = atlas::array::make_view<double, 2>(field);
        if (fsetVec[ji].has(groups_[jg].varInModelFile_)) {
          const atlas::Field inputField = fsetVec[ji][groups_[jg].varInModelFile_];
          const auto inputView = atlas::array::make_view<double, 2>(inputField);
          view.assign(inputView);
        } else {
          throw eckit::Exception("missing " + groups_[jg].varInModelFile_ + " in model file",
            Here());
        }
        rv_->add(field);
      }
    }
  }

  oops::Log::trace() << classname() << "::setReadFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::directCalibration(const oops::FieldSets &) {
  oops::Log::trace() << classname() << "::calibration starting" << std::endl;

  // Get number of layers
  ASSERT(params_.nLayers.value() != boost::none);
  size_t nLayers = *params_.nLayers.value();

  // Allocate data
  data_.resize(groups_.size());
  weight_.resize(nLayers);
  normalization_.resize(nLayers);
  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    // Create layers
    std::vector<std::unique_ptr<LayerBase>> layers;
    for (size_t jBin = 0; jBin < nLayers; ++jBin) {
      data_[jg].emplace_back(LayerFactory::create(params_, fieldsMetaData_, gdata_,
        groups_[jg].name_, groups_[jg].variables_, nx0_, ny0_, groups_[jg].nz0_));
    }
  }
  for (size_t jBin = 0; jBin < nLayers; ++jBin) {
    normalization_[jBin].reset(new oops::FieldSet3D(validTime_, comm_));
  }

  // Setup length-scales
  setupLengthScales();

  // Setup weight
  setupWeight();

  // Setup vertical coordinate
  setupVerticalCoord();

  // Setup resolution
  setupResolution();

  // Setup reduction factors
  setupReductionFactors();

  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    oops::Log::info() << "Info     : Setup of group " << groups_[jg].name_ << ":" << std::endl;

    // Print bins
    oops::Log::info() << "Info     : - Horizontal bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[jg][jBin]->rh() << " < ";
    }
    oops::Log::info() << data_[jg][weight_.size()-1]->rh() << std::endl;
    oops::Log::info() << "Info     : - Vertical bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[jg][jBin]->rv() << " < ";
    }
    oops::Log::info() << data_[jg][weight_.size()-1]->rv() << std::endl;

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      oops::Log::info() << "Info     : * Bin #" << (jBin+1) << ":" << std::endl;

      // Print normalized vertical coordinate
      oops::Log::info() << "Info     :     Normalized vertical coordinate: ";
      for (size_t jz0 = 0; jz0 < groups_[jg].nz0_-1; ++jz0) {
        oops::Log::info() << data_[jg][jBin]->normVertCoord()[jz0] << " < ";
      }
      oops::Log::info() << data_[jg][jBin]->normVertCoord()[groups_[jg].nz0_-1] << std::endl;

      // Setup interpolation
      data_[jg][jBin]->setupInterpolation();

      // Setup kernels
      data_[jg][jBin]->setupKernels();

      // Setup parallelization
      data_[jg][jBin]->setupParallelization();

      // Setup normalization
      data_[jg][jBin]->setupNormalization();
      normalization_[jBin]->add(data_[jg][jBin]->norm()[groups_[jg].name_]);
    }
  }

  // Get square-root of weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    weight_[jBin]->sqrt();
  }

  oops::Log::trace() << classname() << "::calibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  if (comm_.rank() == 0) {
    ASSERT(params_.dataFile.value() != boost::none);

    // NetCDF ids
    int retval, ncid, grpGrpId, layerGrpId;

    // NetCDF file path
    std::string ncfilepath = *params_.dataFile.value();
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Reading file: " << ncfilepath << std::endl;

    // Open file
    if ((retval = nc_open(ncfilepath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

    // Definition mode
    int nLayers;
    if ((retval = nc_get_att_int(ncid, NC_GLOBAL, "nLayers", &nLayers))) ERR(retval);
    ASSERT(nLayers == static_cast<int>(weight_.size()));

    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      // Get group group
      if ((retval = nc_inq_grp_ncid(ncid, groups_[jg].name_.c_str(), &grpGrpId))) ERR(retval);

      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Get layer group
        std::string layerGrpName = "layer_" + std::to_string(jBin);
        if ((retval = nc_inq_grp_ncid(grpGrpId, layerGrpName.c_str(), &layerGrpId))) ERR(retval);
        data_[jg][jBin]->read(layerGrpId);
      }
    }
  }

  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      // Broadcast data
      data_[jg][jBin]->broadcast();
    }
  }

  // Setup reduction factors
  setupReductionFactors();

  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    oops::Log::info() << "Info     : Setup of group " << groups_[jg].name_ << ":" << std::endl;

    // Print bins
    oops::Log::info() << "Info     : - Horizontal bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[jg][jBin]->rh() << " < ";
    }
    oops::Log::info() << data_[jg][weight_.size()-1]->rh() << std::endl;
    oops::Log::info() << "Info     : - Vertical bins: ";
    for (size_t jBin = 0; jBin < weight_.size()-1; ++jBin) {
      oops::Log::info() << data_[jg][jBin]->rv() << " < ";
    }
    oops::Log::info() << data_[jg][weight_.size()-1]->rv() << std::endl;

    for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
      oops::Log::info() << "Info     : * Bin #" << (jBin+1) << ":" << std::endl;

      // Print normalized vertical coordinate
      oops::Log::info() << "Info     :     Normalized vertical coordinate: ";
      oops::Log::info() << data_[jg][jBin]->normVertCoord()[groups_[jg].nz0_-1] << std::endl;

      // Setup interpolation
      data_[jg][jBin]->setupInterpolation();

      // Setup parallelization
      data_[jg][jBin]->setupParallelization();
    }
  }

  // Get square-root of weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    weight_[jBin]->sqrt();
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (comm_.rank() == 0 && (params_.dataFile.value() != boost::none)) {
    // NetCDF ids
    int retval, ncid, grpGrpId, layerGrpId;
    std::vector<std::array<int, 8>> grpIdsVec;

    // NetCDF file path
    std::string ncfilepath = *params_.dataFile.value();
    ncfilepath.append(".nc");
    oops::Log::info() << "Info     : Writing file: " << ncfilepath << std::endl;

    // Create file
    if ((retval = nc_create(ncfilepath.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid))) ERR(retval);

    // Definition mode
    int nLayers = weight_.size();
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "nLayers", NC_INT, 1, &nLayers))) ERR(retval);

    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      // Create group group
      if ((retval = nc_def_grp(ncid, groups_[jg].name_.c_str(), &grpGrpId))) ERR(retval);

      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create layer group
        std::string layerGrpName = "layer_" + std::to_string(jBin);
        if ((retval = nc_def_grp(grpGrpId, layerGrpName.c_str(), &layerGrpId))) ERR(retval);
        std::array<int, 8> grpIds = data_[jg][jBin]->writeDef(layerGrpId);
        grpIdsVec.push_back(grpIds);
      }
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval);

    // Data mode
    size_t jj = 0;
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create layer group
        data_[jg][jBin]->writeData(grpIdsVec[jj]);
        ++jj;
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> FastLAM::fieldsToWrite() const {
  oops::Log::trace() << classname() << "::fieldsToWrite starting" << std::endl;

  // Create vector of pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> pairs;

  // Get vector of configurations
  std::vector<eckit::LocalConfiguration> outputModelFilesConf
    = params_.outputModelFilesConf.value().get_value_or({});

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  for (const auto & conf : outputModelFilesConf) {
    // Get file configuration
    eckit::LocalConfiguration file = getFileConf(comm_, conf);

    // Get parameter
    const std::string param = conf.getString("parameter");

    if (param == "normalized horizontal length-scale") {
      // Create FieldSet3D
      oops::FieldSet3D fset(validTime_, comm_);

      // Copy fields
      for (const auto & var : activeVars_) {
        // Default: missing value
        const size_t nz0 = var.getLevels();
        atlas::Field field = gdata_.functionSpace().createField<double>(
          atlas::option::name(var.name()) | atlas::option::levels(nz0));
        auto view = atlas::array::make_view<double, 2>(field);
        view.assign(util::missingValue<double>());
        fset.add(field);

        for (size_t jg = 0; jg < groups_.size(); ++jg) {
          // Copy field
          if (groups_[jg].varInModelFile_ == var.name()) {
            const atlas::Field rhField = (*rh_)[groups_[jg].name_];
            const auto rhView = atlas::array::make_view<double, 2>(rhField);
            view.assign(rhView);
          }
        }
      }

      // Add pair
      fset.name() = param;
      pairs.push_back(std::make_pair(file, fset));
    }
    if (param == "weight") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create FieldSet3D
        oops::FieldSet3D fset(validTime_, comm_);

        // Copy fields
        for (const auto & var : activeVars_) {
          // Default: missing value
          const size_t nz0 = var.getLevels();
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var.name()) | atlas::option::levels(nz0));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (size_t jg = 0; jg < groups_.size(); ++jg) {
            // Copy field
            if (groups_[jg].varInModelFile_ == var.name()) {
              const atlas::Field wgtSqrtField = (*weight_[jBin])[groups_[jg].name_];
              const auto wgtSqrtView = atlas::array::make_view<double, 2>(wgtSqrtField);
              for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
                if (ghostView(jnode0) == 0) {
                  for (size_t jz0 = 0; jz0 < nz0; ++jz0) {
                    view(jnode0, jz0) = wgtSqrtView(jnode0, jz0)*wgtSqrtView(jnode0, jz0);
                  }
                }
              }
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        fset.name() = param + " - " + std::to_string(jBin);
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
    if (param == "normalization") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create FieldSet3D
        oops::FieldSet3D fset(validTime_, comm_);

        // Copy fields
        for (const auto & var : activeVars_) {
          // Default: missing value
          const size_t nz0 = var.getLevels();
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var.name()) | atlas::option::levels(nz0));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (size_t jg = 0; jg < groups_.size(); ++jg) {
            // Copy field
            if (groups_[jg].varInModelFile_ == var.name()) {
              const atlas::Field normField = (*normalization_[jBin])[groups_[jg].name_];
              const auto normView = atlas::array::make_view<double, 2>(normField);
              view.assign(normView);
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        fset.name() = param + " - " + std::to_string(jBin);
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
    if (param == "normalization accuracy") {
      for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
        // Create FieldSet3D
        oops::FieldSet3D fset(validTime_, comm_);

        // Copy fields
        for (const auto & var : activeVars_) {
          // Default: missing value
          const size_t nz0 = var.getLevels();
          atlas::Field field = gdata_.functionSpace().createField<double>(
            atlas::option::name(var.name()) | atlas::option::levels(nz0));
          auto view = atlas::array::make_view<double, 2>(field);
          view.assign(util::missingValue<double>());
          fset.add(field);

          for (size_t jg = 0; jg < groups_.size(); ++jg) {
            // Copy field
            if (groups_[jg].varInModelFile_ == var.name()) {
              const atlas::Field normField = data_[jg][jBin]->normAcc()[groups_[jg].name_];
              const auto normView = atlas::array::make_view<double, 2>(normField);
              view.assign(normView);
            }
          }
        }

        // Update configuration
        eckit::LocalConfiguration fileCmp(file);
        util::seekAndReplace(fileCmp, "%component%", jBin, 2);

        // Add pair
        fset.name() = param + " - " + std::to_string(jBin);
        pairs.push_back(std::make_pair(fileCmp, fset));
      }
    }
  }

  oops::Log::trace() << classname() << "::fieldsToWrite done" << std::endl;
  return pairs;
}

// -----------------------------------------------------------------------------

void FastLAM::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

void FastLAM::setupLengthScales() {
  oops::Log::trace() << classname() << "::setupLengthScales starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Get rh and rv from yaml
  if (!rh_) {
    ASSERT(params_.rhFromYaml.value() != boost::none);
    rh_.reset(new oops::FieldSet3D(validTime_, comm_));
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      // Get yaml value/profile
      std::vector<double> profile;
      for (const auto & vParams : *params_.rhFromYaml.value()) {
        const std::string vGroupName = vParams.group.value();
        if (vGroupName == groups_[jg].name_) {
          if (profile.size() > 0) {
            throw eckit::UserError("group" + groups_[jg].name_ + " present more that once", Here());
          }
          profile.resize(groups_[jg].nz0_);
          if (vParams.value.value() != boost::none
            && vParams.profile.value() != boost::none) {
            throw eckit::UserError("both value and profile present in the same item", Here());
          }
          if (vParams.value.value() != boost::none) {
            // Copy value
            std::fill(profile.begin(), profile.end(), *vParams.value.value());
          } else if (vParams.profile.value() != boost::none) {
            // Copy profile
            ASSERT(vParams.profile.value()->size() == groups_[jg].nz0_);
            profile = *vParams.profile.value();
          } else {
            throw eckit::UserError("missing value or profile for group " + groups_[jg].name_,
              Here());
          }
        }
      }
      ASSERT(profile.size() > 0);

      // Copy to rh_
      atlas::Field rhField = gdata_.functionSpace().createField<double>(
        atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
      auto rhView = atlas::array::make_view<double, 2>(rhField);
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (ghostView(jnode0) == 0) {
          for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
            rhView(jnode0, jz0) = profile[jz0];
          }
        }
      }
      rh_->add(rhField);
    }
  }

  if (!rv_) {
    ASSERT(params_.rvFromYaml.value() != boost::none);
    rv_.reset(new oops::FieldSet3D(validTime_, comm_));
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      // Get yaml value/profile
      std::vector<double> profile;
      for (const auto & vParams : *params_.rvFromYaml.value()) {
        const std::string vGroupName = vParams.group.value();
        if (vGroupName == groups_[jg].name_) {
          if (profile.size() > 0) {
            throw eckit::UserError("group" + groups_[jg].name_ + " present more that once", Here());
          }
          profile.resize(groups_[jg].nz0_);
          if (vParams.value.value() != boost::none
            && vParams.profile.value() != boost::none) {
            throw eckit::UserError("both value and profile present in the same item", Here());
          }
          if (vParams.value.value() != boost::none) {
            // Copy value
            std::fill(profile.begin(), profile.end(), *vParams.value.value());
          } else if (vParams.profile.value() != boost::none) {
            // Copy profile
            ASSERT(vParams.profile.value()->size() == groups_[jg].nz0_);
            profile = *vParams.profile.value();
          } else {
            throw eckit::UserError("missing value or profile for group " + groups_[jg].name_,
              Here());
          }
        }
      }
      ASSERT(profile.size() > 0);

      // Copy to rv_
      atlas::Field rvField = gdata_.functionSpace().createField<double>(
        atlas::option::name(groups_[jg].name_) | atlas::option::levels(groups_[jg].nz0_));
      auto rvView = atlas::array::make_view<double, 2>(rvField);
      for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
        if (ghostView(jnode0) == 0) {
          for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
            rvView(jnode0, jz0) = profile[jz0];
          }
        }
      }
      rv_->add(rvField);
    }
  }

  // Compute correlation
  oops::Log::info() << "Info     : Correlation rh / rv:" << std::endl;
  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    // Get fields
    atlas::Field rhField = (*rh_)[groups_[jg].name_];
    atlas::Field rvField = (*rv_)[groups_[jg].name_];
    auto rhView = atlas::array::make_view<double, 2>(rhField);
    auto rvView = atlas::array::make_view<double, 2>(rvField);

    // Compute mean
    double rhMean = 0.0;
    double rvMean = 0.0;
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
          rhMean += rhView(jnode0, jz0);
          rvMean += rvView(jnode0, jz0);
        }
      }
    }
    comm_.allReduceInPlace(rhMean, eckit::mpi::sum());
    comm_.allReduceInPlace(rvMean, eckit::mpi::sum());
    rhMean *= 1.0/static_cast<double>(nx0_*ny0_*groups_[jg].nz0_);
    rvMean *= 1.0/static_cast<double>(nx0_*ny0_*groups_[jg].nz0_);

    // Compute horizontal variance and covariance
    double rhVar = 0.0;
    double rvVar = 0.0;
    double cov = 0.0;
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
          rhVar += (rhView(jnode0, jz0)-rhMean)*(rhView(jnode0, jz0)-rhMean);
          rvVar += (rvView(jnode0, jz0)-rvMean)*(rvView(jnode0, jz0)-rvMean);
          cov += (rhView(jnode0, jz0)-rhMean)*(rvView(jnode0, jz0)-rvMean);
        }
      }
    }
    comm_.allReduceInPlace(rhVar, eckit::mpi::sum());
    comm_.allReduceInPlace(rvVar, eckit::mpi::sum());
    comm_.allReduceInPlace(cov, eckit::mpi::sum());

    // Compute horizontal correlation
    double cor = 0.0;
    if (rhVar > 0.0 && rvVar > 0.0) {
      cor = cov/std::sqrt(rhVar*rvVar);
    }

    // Print result
    oops::Log::info() << "Info     : - Group " << groups_[jg].name_ << ": " << cor << std::endl;
  }

  oops::Log::trace() << classname() << "::setupLengthScales done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::setupWeight() {
  oops::Log::trace() << classname() << "::setupWeight starting" << std::endl;

  // Ghost points
  const auto ghostView = atlas::array::make_view<int, 1>(gdata_.functionSpace().ghost());

  // Get function space and grid
  const atlas::functionspace::StructuredColumns fs(gdata_.functionSpace());
  const atlas::StructuredGrid grid(fs.grid());

  // Index fields
  const auto indexX0View = atlas::array::make_view<int, 1>(fs.index_i());
  const auto indexY0View = atlas::array::make_view<int, 1>(fs.index_j());

  // Normalize rh
  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    atlas::Field rhField = (*rh_)[groups_[jg].name_];
    auto rhView = atlas::array::make_view<double, 2>(rhField);

    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        // Cell area
        const int i0 = indexX0View(jnode0)-1;
        const int j0 = indexY0View(jnode0)-1;
        const int im = std::min(std::max(i0-1, 0), static_cast<int>(nx0_-1));
        const int ip = std::min(std::max(i0+1, 0), static_cast<int>(nx0_-1));
        const int jm = std::min(std::max(j0-1, 0), static_cast<int>(ny0_-1));
        const int jp = std::min(std::max(j0+1, 0), static_cast<int>(ny0_-1));
        const double dx = atlas::util::Earth().distance(grid.lonlat(ip, j0), grid.lonlat(im, j0))
          /static_cast<double>(ip-im);
        const double dy = atlas::util::Earth().distance(grid.lonlat(i0, jp), grid.lonlat(i0, jm))
          /static_cast<double>(jp-jm);

        // Normalize rh with cell area square-root
        for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
          rhView(jnode0, jz0) = rhView(jnode0, jz0)/std::sqrt(dx*dy);
        }
      }
    }
  }

  // Initialize weight
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    weight_[jBin].reset(new oops::FieldSet3D(*rh_));
    weight_[jBin]->zero();
  }

  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    // Get rh field general properties
    atlas::Field rhField = (*rh_)[groups_[jg].name_];
    auto rhView = atlas::array::make_view<double, 2>(rhField);
    double minRh = std::numeric_limits<double>::max();
    double maxRh = std::numeric_limits<double>::min();
    for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
      if (ghostView(jnode0) == 0) {
        for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
          minRh = std::min(minRh, rhView(jnode0, jz0));
          maxRh = std::max(maxRh, rhView(jnode0, jz0));
        }
      }
    }
    comm_.allReduceInPlace(minRh, eckit::mpi::min());
    comm_.allReduceInPlace(maxRh, eckit::mpi::max());

    if (weight_.size() > 1) {
      // Prepare binning
      if (minRh < maxRh) {
        double binWidth = (maxRh-minRh)/static_cast<double>(weight_.size()-1);
        for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
          data_[jg][jBin]->rh() = minRh+static_cast<double>(jBin)*binWidth;
        }

        // Compute weight
        for (size_t jnode0 = 0; jnode0 < nodes0_; ++jnode0) {
          if (ghostView(jnode0) == 0) {
            for (size_t jz0 = 0; jz0 < groups_[jg].nz0_; ++jz0) {
              // Raw weight (difference-based)
              std::vector<double> weight(weight_.size());
              double wgtSum = 0.0;
              for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
                atlas::Field fieldWgt = (*weight_[jBin])[groups_[jg].name_];
                auto wgtView = atlas::array::make_view<double, 2>(fieldWgt);
                const double diff = std::abs(rhView(jnode0, jz0)-data_[jg][jBin]->rh())
                  /(maxRh-minRh);
                wgtView(jnode0, jz0) = std::exp(-4.6*diff);  // Factor 4.6 => minimum weight ~0.01
                wgtSum += wgtView(jnode0, jz0);
              }

              // Normalize weight
              for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
                atlas::Field fieldWgt = (*weight_[jBin])[groups_[jg].name_];
                auto wgtView = atlas::array::make_view<double, 2>(fieldWgt);
                wgtView(jnode0, jz0) /= wgtSum;
              }
            }
          }
        }
      } else {
        // Homogeneous rh, using equal weights for each bin
        for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
          data_[jg][jBin]->rh() = minRh;
          atlas::Field fieldWgt = (*weight_[jBin])[groups_[jg].name_];
          auto wgtView = atlas::array::make_view<double, 2>(fieldWgt);
          wgtView.assign(1.0/static_cast<double>(weight_.size()));
        }
      }
    } else {
      // Single bin
      data_[jg][0]->rh() = 0.5*(minRh+maxRh);

      // Set weight to one
      atlas::Field fieldWgt = (*weight_[0])[groups_[jg].name_];
      auto wgtView = atlas::array::make_view<double, 2>(fieldWgt);
      wgtView.assign(1.0);
    }
  }

  oops::Log::trace() << classname() << "::setupWeight done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::setupVerticalCoord() {
  oops::Log::trace() << classname() << "::setupVerticalCoord starting" << std::endl;

  // Setup vertical coordinate
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      data_[jg][jBin]->setupVerticalCoord((*rv_)[groups_[jg].name_],
        (*weight_[jBin])[groups_[jg].name_]);
    }
  }

  oops::Log::trace() << classname() << "::setupVerticalCoord done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::setupResolution() {
  oops::Log::trace() << classname() << "::setupResolution starting" << std::endl;

  // Copy resolution
  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      ASSERT(params_.resol.value() != boost::none);
      data_[jg][jBin]->resol() = *params_.resol.value();
    }
  }

  oops::Log::trace() << classname() << "::setupResolution done" << std::endl;
}

// -----------------------------------------------------------------------------

void FastLAM::setupReductionFactors() {
  oops::Log::trace() << classname() << "::setupReductionFactors starting" << std::endl;

  for (size_t jBin = 0; jBin < weight_.size(); ++jBin) {
    // Define reduction factors
    for (size_t jg = 0; jg < groups_.size(); ++jg) {
      data_[jg][jBin]->rfh() = std::max(data_[jg][jBin]->rh()/data_[jg][jBin]->resol(), 1.0);
      data_[jg][jBin]->rfv() = std::max(data_[jg][jBin]->rv()/data_[jg][jBin]->resol(), 1.0);
    }

    if (params_.strategy.value() == "crossed") {
      // Crossed multivariate strategy: same control variable resolution for all groups
      double rfhMin = 1.0;
      double rfvMin = 1.0;
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        rfhMin = std::min(data_[jg][jBin]->rfh(), rfhMin);
        rfvMin = std::min(data_[jg][jBin]->rfv(), rfvMin);
      }
      for (size_t jg = 0; jg < groups_.size(); ++jg) {
        data_[jg][jBin]->rfh() = rfhMin;
        data_[jg][jBin]->rfv() = rfvMin;
      }
    }
  }

  oops::Log::trace() << classname() << "::setupReductionFactors done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t FastLAM::getGroupIndex(const std::string & var) const {
  oops::Log::trace() << classname() << "::getGroupIndex starting" << std::endl;

  // Initialization
  size_t groupIndex = 0;

  for (size_t jg = 0; jg < groups_.size(); ++jg) {
    // Check group variables
    if (std::find(groups_[jg].variables_.begin(), groups_[jg].variables_.end(), var)
      != groups_[jg].variables_.end()) {
      break;
    }

    // Update group index
    ++groupIndex;
  }

  // Check group index has been found
  ASSERT(groupIndex < groups_.size());

  oops::Log::trace() << classname() << "::getGroupIndex done" << std::endl;
  return groupIndex;
}

// -----------------------------------------------------------------------------

size_t FastLAM::getZ0Offset(const std::string & var) const {
  oops::Log::trace() << classname() << "::getZ0Offset starting" << std::endl;

  // Default value
  size_t z0Offset = 0;

  if (active2dVars_.has(var)) {
    // This is a 2d variable
    if (params_.lev2d.value() == "last") {
      // Use the last level of the group
      z0Offset = groups_[getGroupIndex(var)].nz0_-1;
    }
  }

  oops::Log::trace() << classname() << "::getZ0Offset starting" << std::endl;
  return z0Offset;
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration FastLAM::getFileConf(const eckit::mpi::Comm & comm,
                                               const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::getFileConf starting" << std::endl;

  // Get number of MPI tasks and OpenMP threads
  std::string mpi(std::to_string(comm.size()));
  std::string omp("1");
#ifdef _OPENMP
  # pragma omp parallel
  {
    omp = std::to_string(omp_get_num_threads());
  }
#endif

  // Get IO configuration
  eckit::LocalConfiguration file = conf.getSubConfiguration("file");

  // Replace patterns
  util::seekAndReplace(file, "_MPI_", mpi);
  util::seekAndReplace(file, "_OMP_", omp);

  oops::Log::trace() << classname() << "::getFileConf done" << std::endl;
  return file;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
