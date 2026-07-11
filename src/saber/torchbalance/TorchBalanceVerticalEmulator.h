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

/// Setup vertical Torch emulator: full-column input profiles, compact vertical Jacobian output.
/// Level counts are determined from Atlas field shapes in xb - no input_levels/output_levels
/// TorchScript attributes required.
/// jac_physical must accept (inputs, mask, requested_row_indices, requested_col_indices) and
/// return [nnodes, nRequestedPairs], one value per requested row/column pair.
/// \param xb Trajectory
/// \param comm MPI communicator
/// \param torchscript_path Path to TorchScript model file
/// \param jacFieldSet Output Jacobian fields, each allocated with compact vertical levels:
///        nInputLevels for single-level outputs, otherwise min(nOutputLevels, nInputLevels).
/// \param maskVariable Optional masking variable name (empty string if no masking)
/// \param maskLevel Optional masking variable level (ignored if maskVariable is empty)
/// \return Map of Jacobian field names to (nInputLevels, nOutputLevels) pairs
std::map<std::string, std::pair<int, int>> setupVerticalEmulator(
    const oops::FieldSet3D & xb,
    const eckit::mpi::Comm & comm,
    const std::string & torchscript_path,
    atlas::FieldSet & jacFieldSet,
    const std::string & maskVariable = "",
    int maskLevel = 0);

}  // namespace saber
