/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/GeometryData.h"

#include "saber/fastlam/FastLAMParametersBase.h"
#include "saber/fastlam/Layer.h"

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class Layers {
 public:
  static const std::string classname() {return "saber::fastlam::Layers";}

  typedef FastLAMParametersBase ParametersBase_;

  // Constructors
  Layers() {}
  Layers(const ParametersBase_ & params,
         const oops::GeometryData &,
         const std::string &,
         const size_t &,
         const size_t &,
         const size_t &,
         const size_t &);

  // Accessors
  Layer & operator[](const int ii) {return layers_[ii];}
  const Layer & operator[](const int ii) const {return layers_[ii];}

 private:
  std::vector<Layer> layers_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
