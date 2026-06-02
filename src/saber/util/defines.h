/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace saber {

#ifdef ATLAS_VERSION_46_OR_GREATER
constexpr bool SABER_ATLAS_VERSION_46_OR_GREATER = true;
#else
constexpr bool SABER_ATLAS_VERSION_46_OR_GREATER = false;
#endif  // ifdef ATLAS_VERSION_46_OR_GREATER

}  // namespace saber

