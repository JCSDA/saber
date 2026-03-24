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

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/Diffusion.h"

#include "saber/diffusion/DiffusionParameters.h"

namespace saber {
namespace diffusion {

// --------------------------------------------------------------------------------------

struct Group {
  oops::Variables vars;
  atlas::FieldSet normalization;
  std::unique_ptr<oops::Diffusion> diffusion;
  bool vtDuplicated = false;
  bool varDuplicated = false;
};

// --------------------------------------------------------------------------------------

void createScales(const oops::GeometryData &,
                  std::queue<atlas::Field> &,
                  const std::string &,
                  const eckit::LocalConfiguration &,
                  atlas::FieldSet &);

// ----------------------------------------------------------------------------------

void applyNormSqrt(const atlas::FunctionSpace &,
                   const Group &,
                   atlas::Field &);

// ----------------------------------------------------------------------------------

void addConf(const oops::OptionalParameter<eckit::LocalConfiguration> &,
             std::vector<std::pair<std::string, eckit::LocalConfiguration>> &);

// --------------------------------------------------------------------------------------

void randomize(const oops::GeometryData &,
               const std::vector<Group> &,
               oops::FieldSet3D &);

// --------------------------------------------------------------------------------------

void multiply(const oops::GeometryData &,
              const std::vector<Group> &,
              oops::FieldSet3D &);

// --------------------------------------------------------------------------------------

void filter(const std::vector<Group> &,
            oops::FieldSet3D &);

// --------------------------------------------------------------------------------------

void read(const oops::GeometryData &,
          std::vector<Group> &,
          const DiffusionParameters &);

// --------------------------------------------------------------------------------------

std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs(
  const DiffusionParameters &);

// --------------------------------------------------------------------------------------

void setReadFields(const std::vector<oops::FieldSet3D> &,
                   std::queue<atlas::Field> &);

// --------------------------------------------------------------------------------------

void directCalibration(const oops::GeometryData &,
                       std::vector<Group> &,
                       std::queue<atlas::Field> &,
                       const DiffusionParameters &);

// --------------------------------------------------------------------------------------

}  // namespace diffusion
}  // namespace saber
