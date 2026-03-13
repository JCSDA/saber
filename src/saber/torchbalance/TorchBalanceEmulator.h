/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
* (C) Copyright 2025- UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <map>
#include <string>
#include <utility>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"

namespace saber {

/// Setup Torch emulator: load TorchScript model, run inference, compute Jacobian
/// \param xb Trajectory
/// \param comm MPI communicator
/// \param torchscript_path Path to TorchScript model file
/// \param jacFieldSet Output Jacobian fields (will be filled with computed Jacobians)
/// \param maskVariable Optional masking variable name (empty string if no masking)
/// \param maskLevel Optional masking variable level (ignored if maskVariable is empty)
/// \return Map of Jacobian field names to (input_level, output_level) pairs
std::map<std::string, std::pair<int, int>> setupEmulator(const oops::FieldSet3D & xb,
                                                         const eckit::mpi::Comm & comm,
                                                         const std::string & torchscript_path,
                                                         atlas::FieldSet & jacFieldSet,
                                                         const std::string & maskVariable = "",
                                                         int maskLevel = 0);

}  // namespace saber
