/*
 * (C) Copyright 2025- UCAR
 * (C) Crown Copyright 2024 Met Office
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oops/base/FieldSet4D.h"
#include "oops/base/FieldSets.h"
#include "oops/base/GeometryData.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/OptionalPolymorphicParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/oops/Utilities.h"

// Forward declarations
namespace oops {
  class FieldSet3D;
}

namespace saber {

// -----------------------------------------------------------------------------

class OffDiagWeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OffDiagWeightParameters, Parameters)

 public:
  oops::RequiredParameter<std::vector<std::string>> varPair{"variables pair", this};
  oops::RequiredParameter<double> weight{"value", this};
};

// -----------------------------------------------------------------------------

class SaberCentralBlockGroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberCentralBlockGroupParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> groupName{"group name", this};
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  oops::RequiredParameter<SaberCentralBlockParametersWrapper>
    centralBlock{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>>
    auxOuterBlocksParams{"auxiliary outer blocks", this};

  // Optional parameter specific to "duplicated"
  oops::Parameter<std::vector<std::string>>
    varsOnExtraLevels{"variables on extra levels", {}, this};

  // Optional parameters specific to "duplicated and weighted" strategy
  oops::Parameter<double> defOffDiagWeight{"default off-diagonal weight", 0.0, this};
  oops::OptionalParameter<std::vector<OffDiagWeightParameters>>
    offDiagWeights{"specific off-diagonal weights", this};

  // Direct access to central block parameters
  const SaberBlockParametersBase & centralBlockParams() const
    {return this->centralBlock.value().blockParams();}
};

// -----------------------------------------------------------------------------

class SaberCentralBlockParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberCentralBlockParameters, Parameters)

 public:
  oops::Parameter<std::string> strategy{"multivariate strategy", "deprecated", this};
  // Single block:
  oops::OptionalPolymorphicParameter<SaberBlockParametersBase, SaberCentralBlockFactory>
    singleBlock{"saber block name", this};
  // Or multiple blocks:
  oops::OptionalParameter<std::vector<SaberCentralBlockGroupParameters>>
    groups{"groups", this};

  // Type of setup
  bool doCalibration() const;
  bool doRead() const;

  // Get active variables
  oops::Variables getActiveVars(const oops::Variables &) const;
};

// -----------------------------------------------------------------------------

class SaberCentralBlock : public util::Printable {
 public:
  SaberCentralBlock(const oops::GeometryData & geometryData,
                    const oops::Variables & outerVars,
                    const eckit::Configuration & covarConf,
                    const SaberCentralBlockParameters & params,
                    const oops::FieldSet3D & xb,
                    const oops::FieldSet3D & fg);
  ~SaberCentralBlock() = default;

  // Application methods

  // Block randomization
  void randomize(oops::FieldSet3D &) const;

  // Block multiplication
  void multiply(oops::FieldSet3D &) const;

  // Setup / calibration methods

  // Read block data
  void read();

  // Direct calibration
  void directCalibration(const oops::FieldSets &);

  // Iterative calibration
  void iterativeCalibrationInit();
  void iterativeCalibrationUpdate(const oops::FieldSet3D &);
  void iterativeCalibrationFinal();

  // Write block data
  void write() const;

  // Square-root formulation
  size_t ctlVecSize() const;
  void randomCtlVec(atlas::Field &, const size_t &) const;
  void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const;
  void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const;

  // Return date/time
  const util::DateTime validTime() const {return validTime_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::Variables &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &) const;

  // Calibrate block
  template <typename MODEL>
  void calibrateBlock(const oops::Geometry<MODEL> &,
                      const oops::Variables &,
                      oops::FieldSet4D &,
                      oops::FieldSet4D &,
                      const eckit::Configuration &,
                      std::shared_ptr<SaberOuterBlockChain>,
                      std::shared_ptr<oops::FieldSets>);

  // Adjoint test
  void adjointTest(const double & tol) const;

  // Square-root test
  void sqrtTest(const double & tol) const;

  bool doCalibration() const {
    return std::any_of(doCalibration_.begin(), doCalibration_.end(), [](bool v)
      { return v; });
  }
  bool doRead() const {
    return std::any_of(doRead_.begin(), doRead_.end(), [](bool v) { return v; });
  }
  bool forceWrite() const {
    return std::any_of(forceWrite_.begin(), forceWrite_.end(), [](bool v) { return v; });}

 private:
  // Geometry data
  const oops::GeometryData & geometryData_;
  // Outer variables
  const oops::Variables outerVars_;
  // Valid time
  const util::DateTime validTime_;
  // Parameters
  const SaberCentralBlockParameters params_;

  // Multivariate strategy
  std::string strategy_;

  // Number of groups
  size_t ngroup_;
  // Group variables
  std::vector<oops::Variables> groupInputVars_;
  // Group names
  std::vector<std::string> groupNames_;
  // Group inner variables (containing a single variable with the group name and the metadata and
  // levels of the reference variable)
  std::vector<oops::Variables> groupInnerVars_;
  // Groups need to be calibrated
  std::vector<bool> doCalibration_;
  // Groups need to read MODEL data
  std::vector<bool> doRead_;
  // Groups need to write MODEL data
  std::vector<bool> forceWrite_;
  // Group auxiliary outer block chain
  std::vector<std::unique_ptr<SaberOuterBlockChain>> groupAuxOuterBlockChains_;
  // Group central block
  std::vector<std::unique_ptr<SaberCentralBlockBase>> groupCentralBlocks_;

  // Weights for the "duplicated and weighted" strategy
  std::vector<Eigen::MatrixXd> wgtSqrt_;
  // First level (for 3D and 2D fields summation)
  std::unordered_map<std::string, int> firstLevel_;

  void print(std::ostream &) const {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralBlock::read(const oops::Geometry<MODEL> & geom,
                             const oops::Variables & vars) {
  oops::Log::trace() << "SaberCentralBlock::read starting" << std::endl;
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    groupCentralBlocks_[igroup]->read(geom, vars);
  }
  oops::Log::trace() << "SaberCentralBlock::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralBlock::write(const oops::Geometry<MODEL> & geom) const {
  oops::Log::trace() << "SaberCentralBlock::write starting" << std::endl;
  for (size_t igroup = 0; igroup < ngroup_; ++igroup) {
    groupCentralBlocks_[igroup]->write(geom);
  }
  oops::Log::trace() << "SaberCentralBlock::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SaberCentralBlock::calibrateBlock(const oops::Geometry<MODEL> & geom,
                                       const oops::Variables & outerVars,
                                       oops::FieldSet4D & fset4dXb,
                                       oops::FieldSet4D & fset4dFg,
                                       const eckit::Configuration & conf,
                                       std::shared_ptr<SaberOuterBlockChain> outerBlockChain,
                                       std::shared_ptr<oops::FieldSets> fsetEns) {
  oops::Log::trace() << "SaberCentralBlock::calibrateBlock starting" << std::endl;

  // Iterative ensemble loading flag
  const bool iterativeEnsembleLoading = conf.getBool("iterative ensemble loading");

  // Block calibration
  if (iterativeEnsembleLoading) {
    // Iterative calibration
    oops::Log::info() << "Info     : Iterative calibration" << std::endl;

    // Initialization
    this->iterativeCalibrationInit();

    // Get ensemble size
    size_t nens = getNensFromConfig(conf);

    for (size_t ie = 0; ie < nens; ++ie) {
      // Read ensemble member
      oops::FieldSet3D fset(fset4dXb[0].validTime(), geom.getComm());
      readEnsembleMember(geom, outerVars, conf, ie, fset);

      // Apply outer blocks inverse (all of them)
      oops::Log::info() << "Info     : Apply outer blocks inverse (all of them)" << std::endl;
      if (outerBlockChain) outerBlockChain->leftInverseMultiply(fset);

      // Use FieldSet in the central block
      oops::Log::info() << "Info     : Use FieldSet in the central block" << std::endl;
      this->iterativeCalibrationUpdate(fset);
    }

    // Finalization
    oops::Log::info() << "Info     : Finalization" << std::endl;
    this->iterativeCalibrationFinal();
  } else {
    // Direct calibration
    oops::Log::info() << "Info     : Direct calibration" << std::endl;
    this->directCalibration(*fsetEns);
  }

  oops::Log::trace() << "SaberCentralBlock::calibrateBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
