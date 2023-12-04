/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/fastlam/Layers.h"

#include "oops/util/Logger.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

Layers::Layers(const ParametersBase_ & params,
               const oops::GeometryData & gdata,
               const std::string & myVar,
               const size_t & nx0,
               const size_t & ny0,
               const size_t & nz0,
               const size_t & nLayers) {
  oops::Log::trace() << classname() << "::Layers starting" << std::endl;

  // Create layers
  for (size_t jBin = 0; jBin < nLayers; ++jBin) {
    Layer layer(params, gdata, myVar, jBin, nx0, ny0, nz0);
    layers_.push_back(layer);
  }

  oops::Log::trace() << classname() << "::Layers done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
