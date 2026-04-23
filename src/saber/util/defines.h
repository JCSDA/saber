/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace saber {

#ifdef ATLAS_SCOPE_ISSUE_RESOLVED
constexpr bool SABER_ATLAS_SCOPE_ISSUE_RESOLVED = true;
#else
constexpr bool SABER_ATLAS_SCOPE_ISSUE_RESOLVED = false;
#endif  // ifdef ATLAS_SCOPE_ISSUE_RESOLVED

}  // namespace saber

