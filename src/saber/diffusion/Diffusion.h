/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/diffusion/DiffusionParameters.h"

// forward declarations
namespace oops {
  class Diffusion;
}

namespace saber {

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

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  void read() override;
  std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const override;
  void setReadFields(const std::vector<oops::FieldSet3D> &) override;
  void directCalibration(const oops::FieldSets &) override;

 private:
  void print(std::ostream &) const override {}

  const oops::GeometryData & geom_;
  const std::shared_ptr<oops::Diffusion::DerivedGeom> diffusionGeom_;
  Parameters_ params_;
  oops::Variables vars_;
  std::queue<atlas::Field> calibrateReadFields_;

  struct Group {
    oops::Variables vars;
    atlas::FieldSet normalization;
    std::unique_ptr<oops::Diffusion> diffusion;
    bool vtDuplicated = false;
    bool varDuplicated = false;
  };
  std::vector<Group> groups_;
};
}  // namespace saber
