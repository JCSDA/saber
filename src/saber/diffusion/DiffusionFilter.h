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

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/diffusion/DiffusionImpl.h"
#include "saber/diffusion/DiffusionParameters.h"

// forward declarations
namespace oops {
  class Diffusion;
}

namespace saber {

// -----------------------------------------------------------------------------

/// The diffusion based saber outer block, used for filtering. Diffusion (explicit
/// diffusion in this case) is best for small correlation lengths. If you have large
/// lengths, you're better off using BUMP_NICAS.
class DiffusionFilter : public saber::SaberOuterBlockBase {
 public:
  static const std::string classname() { return "saber::DiffusionFilter"; }
  typedef DiffusionParameters Parameters_;

  DiffusionFilter(const oops::GeometryData &,
                  const oops::Variables &,
                  const eckit::Configuration &,
                  const Parameters_ &,
                  const oops::FieldSet3D &,
                  const oops::FieldSet3D &);

  const oops::GeometryData & innerGeometryData() const override
    {return outerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return outerVars_;}

  void multiply(oops::FieldSet3D & fset) const override
    {diffusion::filter(groups_, fset);}
  void multiplyAD(oops::FieldSet3D &) const override
    {throw eckit::Exception("No adjoint for filter outer blocks", Here());}

  void read() override
    {diffusion::read(outerGeometryData_, groups_, params_);}
  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override
    {return diffusion::getReadConfs(params_);}
  void setReadFields(const std::vector<oops::FieldSet3D> & fvec) override
    {return diffusion::setReadFields(fvec, calibrateReadFields_);}
  void directCalibration(const oops::FieldSets &) override
    {return diffusion::directCalibration(outerGeometryData_, groups_, calibrateReadFields_,
      params_);}

 private:
  void print(std::ostream &) const override {}

  Parameters_ params_;
  std::queue<atlas::Field> calibrateReadFields_;

  std::vector<diffusion::Group> groups_;
};

// -----------------------------------------------------------------------------

}  // namespace saber
