/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
* (C) Copyright 2025- UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "saber/torchbalance/TorchBalanceVerticalEmulator.h"

#include <torch/script.h>
#include <torch/torch.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace saber {

std::map<std::string, std::pair<int, int>> setupVerticalEmulator(
    const oops::FieldSet3D & xb,
    const eckit::mpi::Comm & comm,
    const std::string & torchscript_path,
    atlas::FieldSet & jacFieldSet,
    const std::string & maskVariable,
    int maskLevel) {
  // Load TorchScript model
  torch::jit::script::Module module;
  module = torch::jit::load(torchscript_path);
  module.eval();
  oops::Log::info() << "setupVerticalEmulator: loaded " << torchscript_path << std::endl;

  // Read input variable names
  std::vector<std::string> inputNames;
  for (const auto& v : module.attr("input_names").toListRef())
    inputNames.push_back(v.toString()->string());
  auto emulatorVarInputs = oops::Variables(inputNames);
  const int inputSize = emulatorVarInputs.size();

  // Read output variable names - must be exactly one variable (multiple levels allowed)
  std::vector<std::string> outputNames;
  for (const auto& v : module.attr("output_names").toListRef())
    outputNames.push_back(v.toString()->string());
  if (outputNames.size() != 1) {
    throw eckit::Exception(
        "setupVerticalEmulator: expected exactly 1 output variable, got "
        + std::to_string(outputNames.size()), Here());
  }
  if (!xb.has(outputNames[0])) {
    throw eckit::Exception(
        "setupVerticalEmulator: missing required output field: " + outputNames[0], Here());
  }
  const int nOutputLevels = xb[outputNames[0]].shape(1);

  // Prepare input views; determine level count and offset for each input variable
  std::vector<atlas::array::ArrayView<double, 2>> inputViews;
  std::vector<int> nLevelsPerInput;
  std::vector<int> inputLevelOffsets;
  inputViews.reserve(inputSize);
  nLevelsPerInput.reserve(inputSize);
  inputLevelOffsets.reserve(inputSize);

  int nTotalInputLevels = 0;
  for (int i = 0; i < inputSize; ++i) {
    const std::string& varName = emulatorVarInputs[i].name();
    if (!xb.has(varName)) {
      throw eckit::Exception(
          "setupVerticalEmulator: missing required input field: " + varName, Here());
    }
    inputViews.push_back(atlas::array::make_view<double, 2>(xb[varName]));
    inputLevelOffsets.push_back(nTotalInputLevels);
    nLevelsPerInput.push_back(inputViews.back().shape(1));
    nTotalInputLevels += nLevelsPerInput.back();
  }
  oops::Log::info() << "setupVerticalEmulator: nTotalInputLevels=" << nTotalInputLevels
                    << std::endl;

  const int nnodes = inputViews[0].shape(0);

  // Pack full vertical profiles: input tensor shape [nnodes, nTotalInputLevels]
  // All input levels are needed for the emulator forward pass regardless of which
  // Jacobians are requested.
  std::vector<float> inputData;
  inputData.reserve(nnodes * nTotalInputLevels);
  for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
    for (int i = 0; i < inputSize; ++i) {
      for (int k = 0; k < nLevelsPerInput[i]; ++k) {
        inputData.push_back(static_cast<float>(inputViews[i](jnode, k)));
      }
    }
  }
  torch::Tensor inputs = torch::from_blob(inputData.data(),
                                          {nnodes, nTotalInputLevels},
                                          torch::kFloat).clone();

  // Prepare mask tensor [nnodes, 1]
  torch::Tensor mask;
  if (!maskVariable.empty()) {
    auto mask_view = atlas::array::make_view<double, 2>(xb[maskVariable]);
    const int mask_nlevels = mask_view.shape(1);
    if (maskLevel < 0 || maskLevel >= mask_nlevels) {
      throw eckit::Exception(
          "setupVerticalEmulator: mask level " + std::to_string(maskLevel)
          + " out of bounds for '" + maskVariable + "' with "
          + std::to_string(mask_nlevels) + " levels", Here());
    }
    std::vector<float> maskData;
    maskData.reserve(nnodes);
    for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode)
      maskData.push_back(static_cast<float>(mask_view(jnode, maskLevel)));
    mask = torch::from_blob(maskData.data(), {nnodes, 1}, torch::kFloat).clone();
  } else {
    mask = torch::ones({nnodes, 1}, torch::kFloat);
  }

  // Resolve each requested Jacobian field to its model input variable, and build
  // the compact list of row/column pairs to pass to jac_physical.
  std::map<std::string, int> inputNameToIdx;
  for (int i = 0; i < inputSize; ++i)
    inputNameToIdx[emulatorVarInputs[i].name()] = i;

  struct FieldRequest {
    atlas::Field field;
    int nInputLevels;
    int nStoredLevels;
    int reqPairOffset;  // starting pair index within the returned compact Jacobian tensor
  };
  std::vector<FieldRequest> fieldRequests;

  std::vector<int64_t> requestedRowIndices;
  std::vector<int64_t> requestedColIndices;
  for (auto& jacField : jacFieldSet) {
    const std::string& fieldName = jacField.name();
    const size_t pos = fieldName.find("div_d");
    if (pos == std::string::npos) continue;
    const std::string inputSuffix = fieldName.substr(pos + 5);

    auto it = inputNameToIdx.find(inputSuffix);
    if (it == inputNameToIdx.end()) {
      oops::Log::warning() << "setupVerticalEmulator: could not match field "
                           << fieldName << " to any input variable" << std::endl;
      continue;
    }

    const int inputVarIdx  = it->second;
    const int levelOffset  = inputLevelOffsets[inputVarIdx];
    const int nInputLevels = nLevelsPerInput[inputVarIdx];
    const int nStoredLevels = static_cast<int>(jacField.shape(1));

    FieldRequest fr;
    fr.field        = jacField;
    fr.nInputLevels = nInputLevels;
    fr.nStoredLevels = nStoredLevels;
    fr.reqPairOffset = static_cast<int>(requestedColIndices.size());
    fieldRequests.push_back(fr);

    for (int k = 0; k < nStoredLevels; ++k) {
      requestedRowIndices.push_back(nOutputLevels == 1 ? 0 : k);
      requestedColIndices.push_back(levelOffset + k);
    }
  }

  oops::Log::info() << "setupVerticalEmulator: requesting "
                    << requestedColIndices.size() << " compact row/column pairs" << std::endl;

  // Call jac_physical(inputs, mask, requested_row_indices, requested_col_indices).
  // Returns only the requested compact Jacobian entries: [nnodes, nRequestedPairs].
  torch::Tensor reqRowTensor = torch::tensor(requestedRowIndices, torch::kLong);
  torch::Tensor reqColTensor = torch::tensor(requestedColIndices, torch::kLong);
  std::vector<torch::jit::IValue> inputs_vec{inputs, mask, reqRowTensor, reqColTensor};
  torch::Tensor jac = module.get_method("jac_physical")(inputs_vec).toTensor();

  const int nRequestedPairs = static_cast<int>(requestedColIndices.size());
  oops::Log::info() << "setupVerticalEmulator: nOutputLevels=" << nOutputLevels << std::endl;

  if (jac.dim() != 2 || jac.size(1) != nRequestedPairs) {
    throw eckit::Exception(
        "setupVerticalEmulator: jac_physical must return [nnodes, "
        + std::to_string(nRequestedPairs) + "] compact entries", Here());
  }

  auto jac_accessor = jac.accessor<float, 2>();
  std::map<std::string, std::pair<int, int>> levelMetadata;

  // Unpack: each field's compact block starts at fr.reqPairOffset.
  for (auto& fr : fieldRequests) {
    auto jac_view = atlas::array::make_view<double, 2>(fr.field);
    for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
      for (int k = 0; k < fr.nStoredLevels; ++k)
        jac_view(jnode, k) = static_cast<double>(jac_accessor[jnode][fr.reqPairOffset + k]);
    }
    fr.field.haloExchange();

    levelMetadata[fr.field.name()] = std::make_pair(fr.nInputLevels, nOutputLevels);
    oops::Log::info() << "  setupVerticalEmulator: " << fr.field.name()
                      << " nInputLevels=" << fr.nInputLevels
                      << " nOutputLevels=" << nOutputLevels << std::endl;
  }

  return levelMetadata;
}

}  // namespace saber
