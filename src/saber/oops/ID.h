/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace generic {

// -----------------------------------------------------------------------------

class IDParameters : public SaberCentralBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(IDParameters, SaberCentralBlockParametersBase)
};

// -----------------------------------------------------------------------------

class ID : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::generic::ID";}

  typedef IDParameters Parameters_;

  ID(const oops::GeometryData &,
     const std::vector<size_t> &,
     const oops::Variables &,
     const Parameters_ &,
     const atlas::FieldSet &,
     const atlas::FieldSet &,
     const std::vector<atlas::FieldSet> &);

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace generic
}  // namespace saber
