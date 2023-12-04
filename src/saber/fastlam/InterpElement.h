/*
 * (C) Copyright 2023 Meteorlogisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace saber {
namespace fastlam {

// -----------------------------------------------------------------------------

class InterpElement {
 public:
  static const std::string classname() {return "saber::fastlam::InterpElement";}

  // Constructor
  InterpElement(const size_t & index1,
                const std::vector<std::pair<size_t, double>> & operations) :
    interpType_("n"),
    index1_(index1),
    index2_(0),
    operations_(operations) {}
  InterpElement(const std::string & interpType,
                const size_t & index1,
                const size_t & index2,
                const std::vector<std::pair<size_t, double>> & operations) :
    interpType_(interpType),
    index1_(index1),
    index2_(index2),
    operations_(operations) {}

  // Accessors
  const std::string & interpType() const {return interpType_;}
  const size_t & index1() const {return index1_;}
  const size_t & index2() const {return index2_;}
  const std::vector<std::pair<size_t, double>> & operations() const {return operations_;}

 private:
  // Interpolation type ("c": colocated, "x": linear along x, "y": linear along y, "b": bilinear)
  std::string interpType_;

  // Global indices of the most southwestern point
  size_t index1_;
  size_t index2_;

  // Operations (vector of pairs: index in mVec and weight)
  std::vector<std::pair<size_t, double>> operations_;
};

// -----------------------------------------------------------------------------

}  // namespace fastlam
}  // namespace saber
