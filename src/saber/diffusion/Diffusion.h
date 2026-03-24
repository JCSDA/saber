/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/diffusion/DiffusionImpl.h"
#include "saber/diffusion/DiffusionParameters.h"

// forward declarations
namespace oops {
  class Diffusion;
}

namespace saber {

// -----------------------------------------------------------------------------

/// The diffusion based correlation/localization saber central block. Diffusion (explicit
/// diffusion in this case) is best for small correlation lengths. If you have large
/// lengths, you're better off using BUMP_NICAS.
class Diffusion : public saber::SaberCentralBlockBase {
 public:
  static const std::string classname() { return "saber::Diffusion"; }
  typedef DiffusionParameters Parameters_;

  Diffusion(const oops::GeometryData &,
            const oops::Variables &,
            const eckit::Configuration &,
            const Parameters_ &,
            const oops::FieldSet3D &,
            const oops::FieldSet3D &);

  void randomize(oops::FieldSet3D & fset) const override
    {diffusion::randomize(geometryData(), groups_, fset);}
  void multiply(oops::FieldSet3D & fset) const override
    {diffusion::multiply(geometryData(), groups_, fset);}

  void read() override
    {diffusion::read(geometryData(), groups_, params_);}
  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override
    {return diffusion::getReadConfs(params_);}
  void setReadFields(const std::vector<oops::FieldSet3D> & fvec) override
    {return diffusion::setReadFields(fvec, calibrateReadFields_);}
  void directCalibration(const oops::FieldSets &) override
    {return diffusion::directCalibration(geometryData(), groups_, calibrateReadFields_, params_);}

  size_t ctlVecSize() const override
    {return ctlVecSize_;}

 private:
  void print(std::ostream &) const override {}

  size_t ctlVecSize_;
  Parameters_ params_;
  std::queue<atlas::Field> calibrateReadFields_;

  std::vector<diffusion::Group> groups_;
};

// -----------------------------------------------------------------------------

}  // namespace saber
