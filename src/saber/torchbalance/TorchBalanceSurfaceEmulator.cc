/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
* (C) Copyright 2025- UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "saber/torchbalance/TorchBalanceSurfaceEmulator.h"

#include <torch/script.h>
#include <torch/torch.h>

#include <sstream>
#include <utility>
#include <vector>

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace saber {

std::map<std::string, std::pair<int, int>> setupSurfaceEmulator(
    const oops::FieldSet3D & xb,
    const eckit::mpi::Comm & comm,
    const std::string & torchscript_path,
    atlas::FieldSet & jacFieldSet,
    const std::string & maskVariable,
    int maskLevel) {
  // Log PyTorch version information
  oops::Log::info() << "PyTorch C++ API version: " << TORCH_VERSION_MAJOR << "."
                    << TORCH_VERSION_MINOR << "." << TORCH_VERSION_PATCH << std::endl;

  // Load TorchScript model
  // ----------------------
  torch::jit::script::Module module;
  module = torch::jit::load(torchscript_path);
  module.eval();
  oops::Log::info() << "Successfully loaded TorchScript module" << std::endl;

  // Get the input information from the torchscript module
  // -----------------------------------------------------
  // input names
  std::vector<std::string> inputNames;
  auto input_iv = module.attr("input_names");
  for (const auto& v : input_iv.toListRef()) {
    inputNames.push_back(v.toString()->string());
  }
  auto emulatorVarInputs = oops::Variables(inputNames);
  int inputSize = emulatorVarInputs.size();

  // input levels
  std::vector<int> inputLevels;
  auto input_lv = module.attr("input_levels");
  for (const auto& v : input_lv.toListRef()) {
    inputLevels.push_back(v.toInt());
  }

  oops::Log::info() << "Inputs ";
  oops::Log::info() << emulatorVarInputs << std::endl;
  oops::Log::info() << "Level indices: " << inputLevels << std::endl;

  // Get the output information from the torchscript module
  // ------------------------------------------------------
  // output names
  std::vector<std::string> outputNames;
  auto output_iv = module.attr("output_names");
  for (const auto& v : output_iv.toListRef()) {
    outputNames.push_back(v.toString()->string());
  }
  auto emulatorVarOutputs = oops::Variables(outputNames);

  // output levels
  std::vector<int> outputLevels;
  auto output_lv = module.attr("output_levels");
  for (const auto& v : output_lv.toListRef()) {
    outputLevels.push_back(v.toInt());
  }

  oops::Log::info() << "Outputs ";
  oops::Log::info() << emulatorVarOutputs << std::endl;
  oops::Log::info() << "Level indices: " << outputLevels << std::endl;

  // We're only expecting one output for now
  if (emulatorVarOutputs.size() != 1) {
    throw eckit::Exception("setupSurfaceEmulator expected exactly 1 output from emulator", Here());
  }

  // Prepare the input field views
  std::vector<atlas::array::ArrayView<double, 2>> inputViews;
  inputViews.reserve(inputSize);
  for (int i = 0; i < inputSize; ++i) {
    const std::string varName = emulatorVarInputs[i].name();
    if (xb.has(varName)) {
      inputViews.push_back(atlas::array::make_view<double, 2>(xb[varName]));
    } else {
      std::ostringstream err;
      err << "setupSurfaceEmulator missing required input field: " << varName;
      throw eckit::Exception(err.str(), Here());
    }
  }

  // Validate input level indices before data collection
  for (int i = 0; i < inputSize; ++i) {
    const std::string varName = emulatorVarInputs[i].name();
    int nlevels = inputViews[i].shape(1);
    int requested_level = inputLevels[i];

    if (requested_level < 0 || requested_level >= nlevels) {
      std::ostringstream err;
      err << "setupSurfaceEmulator: Input level index out of bounds for variable " << varName
          << ". Requested level " << requested_level << " but variable has "
          << nlevels << " levels.";
      throw eckit::Exception(err.str(), Here());
    }
  }

  // Collect input data for all nodes
  const int nnodes = inputViews[0].shape(0);
  std::vector<float> inputData;
  inputData.reserve(nnodes * inputSize);

  for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
    for (int i = 0; i < inputSize; ++i) {
      inputData.push_back(static_cast<float>(inputViews[i](jnode, inputLevels[i])));
    }
  }

  oops::Log::debug() << "Rank " << comm.rank() << ": Processing " << nnodes
                     << " nodes" << std::endl;

  // Create input tensor
  torch::Tensor inputs = torch::from_blob(inputData.data(),
                                          {nnodes, inputSize},
                                          torch::kFloat).clone();

  // Prepare mask tensor (always required by jac_physical)
  torch::Tensor mask;

  if (!maskVariable.empty()) {
    // Use custom mask from specified variable
    oops::Log::info() << "Rank " << comm.rank() << ": Preparing mask from variable '"
                      << maskVariable << "' at level " << maskLevel << std::endl;

    auto mask_view = atlas::array::make_view<double, 2>(xb[maskVariable]);

    // Validate mask level
    int mask_nlevels = mask_view.shape(1);
    if (maskLevel < 0 || maskLevel >= mask_nlevels) {
      std::ostringstream err;
      err << "setupSurfaceEmulator: Mask level " << maskLevel
          << " out of bounds for variable '" << maskVariable
          << "' with " << mask_nlevels << " levels";
      throw eckit::Exception(err.str(), Here());
    }

    // Collect mask data for all nodes
    std::vector<float> maskData;
    maskData.reserve(nnodes);
    for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
      maskData.push_back(static_cast<float>(mask_view(jnode, maskLevel)));
    }

    // Create mask tensor with shape [nnodes, 1]
    mask = torch::from_blob(maskData.data(), {nnodes, 1}, torch::kFloat).clone();
  } else {
    // Default mask: all ones (no masking)
    oops::Log::info() << "Rank " << comm.rank() << ": Using default mask (all ones)" << std::endl;
    mask = torch::ones({nnodes, 1}, torch::kFloat);
  }

  // Call the jac_physical method directly from the model
  // The model's jac_physical() method returns d(output_physical)/d(input_physical)
  // with shape [nnodes, outputSize, inputSize]
  // Mask is always passed as the second argument
  // TODO(Guillaume): Pass the name of the jacobian method through the configuration; mask can
  //       be an ancillary tensor from which a mask can be calculated if needed (rename mask
  //       to ancTensor or similar)
  oops::Log::debug() << "Rank " << comm.rank() <<
    ": Computing Jacobian via model.jac_physical() method..." << std::endl;
  std::vector<torch::jit::IValue> inputs_vec{inputs, mask};
  torch::Tensor jac = module.get_method("jac_physical")(inputs_vec).toTensor();

  // Extract Jacobian values for all fields in jacFieldSet
  // The jac_physical() method returns shape [nnodes, outputSize, inputSize]
  // Already in physical units (no scaling needed)

  oops::Log::debug() << "Extracting Jacobians for "
                     << jacFieldSet.size() << " fields" << std::endl;

  // Get accessor once for efficiency
  auto jac_accessor = jac.accessor<float, 3>();

  // Build level metadata map
  std::map<std::string, std::pair<int, int>> levelMetadata;

  // Process each Jacobian field
  for (auto& jacField : jacFieldSet) {
    const std::string& fieldName = jacField.name();
    oops::Log::info() << "  Processing Jacobian field: " << fieldName << std::endl;

    // Parse the field name to determine which Jacobian this is
    // Expected format: "doutput_div_dinput" (e.g., "dtair_div_dsst")
    // Extract the input variable name after "div_d"

    int input_var_index = -1;
    size_t pos = fieldName.find("div_d");

    if (pos != std::string::npos) {
      // Extract the part after "div_d"
      std::string inputSuffix = fieldName.substr(pos + 5);  // 5 = length of "div_d"

      // Try to match with input variable names
      for (int i = 0; i < inputSize; ++i) {
        std::string inputVarName = emulatorVarInputs[i].name();
        if (inputVarName == inputSuffix) {
          input_var_index = i;
          break;
        }
      }
    }

    oops::Log::debug() << "    Jacobian w.r.t. input variable: "
                      << emulatorVarInputs[input_var_index].name()
                      << " (index " << input_var_index << ")" << std::endl;

    // Store level metadata for this Jacobian field
    // input level from inputLevels[input_var_index]
    // output level from outputLevels[0] (assuming single output)
    levelMetadata[fieldName] = std::make_pair(inputLevels[input_var_index], outputLevels[0]);

    oops::Log::debug() << "    Level metadata: input_level=" << inputLevels[input_var_index]
                      << ", output_level=" << outputLevels[0] << std::endl;

    // Fill in the Jacobian Field
    auto jac_view = atlas::array::make_view<double, 2>(jacField);

    // Fill all nodes with Jacobian values from the emulator
    // The emulator handles masking internally (returns 0 where appropriate)
    for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
      jac_view(jnode, 0) = static_cast<double>(jac_accessor[jnode][0][input_var_index]);
    }

    // Synchronize halo/ghost nodes across MPI ranks
    jacField.haloExchange();
  }

  oops::Log::info() << "All Jacobian fields complete" << std::endl;

  return levelMetadata;
}

}  // namespace saber
