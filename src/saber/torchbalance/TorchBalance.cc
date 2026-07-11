/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
* (C) Copyright 2025- UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <algorithm>

#include "saber/torchbalance/TorchBalance.h"
#include "saber/torchbalance/TorchBalanceSurfaceEmulator.h"
#include "saber/torchbalance/TorchBalanceVerticalEmulator.h"

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"

namespace saber {

// --------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<TorchBalance> makerTorchBalance_("TorchBalance");

// --------------------------------------------------------------------------------------

TorchBalance::TorchBalance(
    const oops::GeometryData & outerGeometryData,
    const oops::Variables & outerVars,
    const eckit::Configuration & tbConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerVars_(outerVars),
    innerGeometryData_(outerGeometryData),
    jac_(nullptr)
{
  oops::Log::trace() << "TorchBalance constructed, inner vars: " << innerVars_ << std::endl;

  // ---- Phase 1: collect all Jacobian field specs (name + level count) ----

  std::vector<std::string> jacStr;
  std::vector<int> jacLevels;

  // Surface emulators: single-level Jacobian per (output, input) pair
  const auto& emulators = params.surfaceEmulators.value();
  for (const auto& emulator : emulators) {
    const std::string& emulatorName = emulator.name.value();
    for (const auto& inputVar : emulator.jacobianWrt.value()) {
      jacStr.push_back("d" + emulatorName + "_div_d" + inputVar);
      jacLevels.push_back(1);
    }
  }

  // Vertical emulators: compact vertical-block Jacobian.
  if (params.verticalEmulators.value() != boost::none) {
    for (const auto& emulator : *params.verticalEmulators.value()) {
      const std::string& emulatorName = emulator.name.value();
      if (!xb.has(emulatorName)) {
        throw eckit::Exception(
            "TorchBalance: vertical emulator output '" + emulatorName
            + "' not found in background state", Here());
      }
      const int nOutputLevels = xb[emulatorName].shape(1);
      for (const auto& inputVar : emulator.jacobianWrt.value()) {
        if (!xb.has(inputVar)) {
          throw eckit::Exception(
              "TorchBalance: vertical emulator input '" + inputVar
              + "' not found in background state", Here());
        }
        const int nInputLevels = xb[inputVar].shape(1);
        jacStr.push_back("d" + emulatorName + "_div_d" + inputVar);
        jacLevels.push_back(nOutputLevels == 1 ? nInputLevels :
                            std::min(nOutputLevels, nInputLevels));
      }
    }
  }

  oops::Log::info() << "TorchBalance: allocating " << jacStr.size()
                    << " Jacobian fields" << std::endl;

  // ---- Phase 2: allocate the Jacobian FieldSet ----

  oops::Variables jacVars(jacStr);
  for (size_t i = 0; i < jacStr.size(); ++i)
    jacVars[jacStr[i]].setLevels(jacLevels[i]);
  const double INITIAL_VALUE = 0.0;
  auto fs = innerGeometryData_.functionSpace();
  jac_ = util::createFieldSet(fs, jacVars, INITIAL_VALUE);

  // ---- Phase 3: run surface emulators ----

  for (const auto& emulator : emulators) {
    const std::string& emulatorName = emulator.name.value();
    const std::string& emulatorPath = emulator.path.value();
    const auto& jacWrt = emulator.jacobianWrt.value();

    oops::Log::info() << "TorchBalance surface emulator: " << emulatorName
                      << "  path: " << emulatorPath << std::endl;

    std::string maskVariable = "";
    int maskLevel = 0;
    if (emulator.maskVariable.value() != boost::none)
      maskVariable = *emulator.maskVariable.value();
    if (emulator.maskLevel.value() != boost::none)
      maskLevel = *emulator.maskLevel.value();

    atlas::FieldSet emulatorJacFieldSet;
    for (const auto& inputVar : jacWrt)
      emulatorJacFieldSet.add(jac_["d" + emulatorName + "_div_d" + inputVar]);

    auto emulatorLevelMetadata = setupSurfaceEmulator(xb, outerGeometryData.comm(),
                              emulatorPath, emulatorJacFieldSet,
                              maskVariable, maskLevel);

    // Surface: denominator_level and numerator_level from emulator;
    //          nInputLevels = nOutputLevels = 1
    for (const auto& entry : emulatorLevelMetadata) {
      jacLevelMetadata_[entry.first] = {entry.second.first, entry.second.second, 1, 1};
      oops::Log::info() << "  surface Jacobian " << entry.first
                        << ": denom_level=" << entry.second.first
                        << " numer_level=" << entry.second.second << std::endl;
    }
  }

  // ---- Phase 4: run vertical emulators ----

  if (params.verticalEmulators.value() != boost::none) {
    for (const auto& emulator : *params.verticalEmulators.value()) {
      const std::string& emulatorName = emulator.name.value();
      const std::string& emulatorPath = emulator.path.value();
      const auto& jacWrt = emulator.jacobianWrt.value();

      oops::Log::info() << "TorchBalance vertical emulator: " << emulatorName
                        << "  path: " << emulatorPath << std::endl;

      std::string maskVariable = "";
      int maskLevel = 0;
      if (emulator.maskVariable.value() != boost::none)
        maskVariable = *emulator.maskVariable.value();
      if (emulator.maskLevel.value() != boost::none)
        maskLevel = *emulator.maskLevel.value();

      atlas::FieldSet emulatorJacFieldSet;
      for (const auto& inputVar : jacWrt)
        emulatorJacFieldSet.add(jac_["d" + emulatorName + "_div_d" + inputVar]);

      auto emulatorLevelMetadata = setupVerticalEmulator(xb, outerGeometryData.comm(),
                                                         emulatorPath, emulatorJacFieldSet,
                                                         maskVariable, maskLevel);

      // Vertical: denom_level = numer_level = 0; nInputLevels and nOutputLevels from metadata
      for (const auto& entry : emulatorLevelMetadata) {
        jacLevelMetadata_[entry.first] = {0, 0, entry.second.first, entry.second.second};
        oops::Log::info() << "  vertical Jacobian " << entry.first
                          << ": nInputLevels=" << entry.second.first
                          << " nOutputLevels=" << entry.second.second << std::endl;
      }
    }
  }

  // Save the Jacobian fields to file if enabled
  if (params.saveJacobians.value()) {
    oops::Log::debug() << "Saving Jacobian FieldSet for debugging" << std::endl;
    // Get grid identifier
    const std::string gridUid = util::getGridUid(outerGeometryData.functionSpace());

    // Create a minimal configuration for debugging writes
    eckit::LocalConfiguration debugConf;
    const std::string jacPath = "jacobian_" + gridUid;
    debugConf.set("filepath", jacPath);

    // Wrap FieldSetImpl pointer in an atlas::FieldSet handle for writing
    atlas::FieldSet jacHandle(jac_);
    util::writeFieldSet(outerGeometryData.comm(), debugConf, jacHandle);
    oops::Log::info() << "Saved Jacobian FieldSet (debug write) to: " << jacPath << std::endl;
  }
}

// --------------------------------------------------------------------------------------

void TorchBalance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << "TorchBalance::multiply start" << std::endl;
  auto ghostView = atlas::array::make_view<int, 1>(innerGeometryData_.functionSpace().ghost());

  // Make a copy of the input fields to read from (avoid reading modified values)
  atlas::FieldSet fsetCopy;
  for (const auto& field : fset.fieldSet()) {
    fsetCopy.add(field.clone());
  }

  // Loop over each variable in the field set
  for (const auto& var : innerVars_) {
    const std::string& varName = var.name();
    if (!fset.has(varName)) continue;

    auto var_view = atlas::array::make_view<double, 2>(fset[varName]);

    // Forward: d(varName) += sum_j { d{varName}_div_d{otherVar} * d{otherVar} }
    // Loop over all Jacobian fields and find those where varName is the numerator
    for (const auto& jacField : jac_) {
      const std::string& jacName = jacField.name();

      // Parse: d{numerator}_div_d{denominator}
      std::string numerator, denominator;
      if (!parseJacobianName(jacName, numerator, denominator)) continue;

      // Only process Jacobians where this variable is the numerator
      if (numerator != varName) continue;
      if (!fset.has(denominator)) continue;

      // Look up level metadata
      int denominator_level, numerator_level, nInputLevels, nOutputLevels;
      if (!getJacobianLevels(jacName, denominator_level, numerator_level,
                             nInputLevels, nOutputLevels)) continue;

      auto denominator_view = atlas::array::make_view<double, 2>(fsetCopy[denominator]);
      auto jac_view         = atlas::array::make_view<double, 2>(jac_[jacName]);

      oops::Log::debug() << "  multiplyFWD: d" << varName << " += "
                         << jacName << " * d" << denominator << std::endl;
      const int nnodes = jac_[jacName].shape(0);

      if (nOutputLevels == 1) {
        for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
          if (ghostView(jnode) > 0) continue;
          for (int k = 0; k < nInputLevels; ++k)
            var_view(jnode, numerator_level) +=
              jac_view(jnode, k) * denominator_view(jnode, denominator_level + k);
        }
      } else {
        const int nDiagonalLevels = std::min(nOutputLevels, nInputLevels);
        for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
          if (ghostView(jnode) > 0) continue;
          for (int k = 0; k < nDiagonalLevels; ++k)
            var_view(jnode, numerator_level + k) +=
              jac_view(jnode, k) * denominator_view(jnode, denominator_level + k);
        }
      }
    }
  }

  oops::Log::trace() << "TorchBalance::multiply done" << std::endl;
}

// --------------------------------------------------------------------------------------

void TorchBalance::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << "TorchBalance::multiplyAD start" << std::endl;
  auto ghostView = atlas::array::make_view<int, 1>(innerGeometryData_.functionSpace().ghost());

  // Make a copy of the input fields to read from (avoid reading modified values)
  atlas::FieldSet fsetCopy;
  for (const auto& field : fset.fieldSet()) {
    fsetCopy.add(field.clone());
  }

  // Loop over each variable in the field set
  for (const auto& var : innerVars_) {
    const std::string& varName = var.name();
    if (!fset.has(varName)) continue;

    auto var_view = atlas::array::make_view<double, 2>(fset[varName]);

    // Adjoint: d(varName) += sum_j { d{otherVar}_div_d{varName} * d{otherVar} }
    // Loop over all Jacobian fields and find those where varName is the denominator
    for (const auto& jacField : jac_) {
      const std::string& jacName = jacField.name();

      // Parse: d{numerator}_div_d{denominator}
      std::string numerator, denominator;
      if (!parseJacobianName(jacName, numerator, denominator)) continue;

      // Only process Jacobians where this variable is the denominator
      if (denominator != varName) continue;
      if (!fset.has(numerator)) continue;

      // Look up level metadata
      int denominator_level, numerator_level, nInputLevels, nOutputLevels;
      if (!getJacobianLevels(jacName, denominator_level, numerator_level,
                             nInputLevels, nOutputLevels)) continue;

      auto numerator_view = atlas::array::make_view<double, 2>(fsetCopy[numerator]);
      auto jac_view       = atlas::array::make_view<double, 2>(jac_[jacName]);

      oops::Log::debug() << "  multiplyAD: d" << varName << " += "
                         << jacName << " * d" << numerator << std::endl;

      const int nnodes = jac_[jacName].shape(0);
      if (nOutputLevels == 1) {
        for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
          if (ghostView(jnode) > 0) continue;
          for (int k = 0; k < nInputLevels; ++k)
            var_view(jnode, denominator_level + k) +=
              jac_view(jnode, k) * numerator_view(jnode, numerator_level);
        }
      } else {
        const int nDiagonalLevels = std::min(nOutputLevels, nInputLevels);
        for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
          if (ghostView(jnode) > 0) continue;
          for (int k = 0; k < nDiagonalLevels; ++k)
            var_view(jnode, denominator_level + k) +=
              jac_view(jnode, k) * numerator_view(jnode, numerator_level + k);
        }
      }
    }
  }

  oops::Log::trace() << "TorchBalance::multiplyAD done" << std::endl;
}

// --------------------------------------------------------------------------------------

void TorchBalance::print(std::ostream &) const {}

// --------------------------------------------------------------------------------------

bool TorchBalance::parseJacobianName(const std::string & jacName,
                                     std::string & numerator,
                                     std::string & denominator) {
  // Expected format: "d{numerator}_div_d{denominator}"
  if (jacName.empty() || jacName[0] != 'd') return false;
  const size_t div_pos = jacName.rfind("_div_d");
  if (div_pos == std::string::npos) return false;
  numerator   = jacName.substr(1, div_pos - 1);
  denominator = jacName.substr(div_pos + 6);
  return !numerator.empty() && !denominator.empty();
}

// --------------------------------------------------------------------------------------

bool TorchBalance::getJacobianLevels(const std::string & jacName,
                                     int & denominator_level,
                                     int & numerator_level,
                                     int & nInputLevels,
                                     int & nOutputLevels) const {
  auto it = jacLevelMetadata_.find(jacName);
  if (it == jacLevelMetadata_.end()) return false;
  denominator_level = it->second.denominator_level;
  numerator_level   = it->second.numerator_level;
  nInputLevels      = it->second.nInputLevels;
  nOutputLevels     = it->second.nOutputLevels;
  return true;
}

// --------------------------------------------------------------------------------------

}  // namespace saber
