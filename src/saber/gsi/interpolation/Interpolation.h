/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/gsi/interpolation/unstructured_interp/UnstructuredInterpolation.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class InterpolationParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, SaberBlockParametersBase)

 public:
  // File containing grid and coefficients
  oops::RequiredParameter<std::string> GSIFile{"gsi error covariance file", this};
  oops::RequiredParameter<std::string> GSINML{"gsi berror namelist file", this};
  oops::RequiredParameter<std::string> GSIVGRD{"gsi akbk", this};

  // Handle vertical top-2-bottom and vice-verse wrt to GSI
  oops::Parameter<bool> vflip{"flip vertical grid", true, this};

  // Processor layout
  oops::Parameter<size_t> layoutx{"processor layout x direction", 1, this};
  oops::Parameter<size_t> layouty{"processor layout y direction", 1, this};

  // Debugging mode
  oops::Parameter<bool> debugMode{"debugging mode", false, this};
  oops::Parameter<bool> bypassGSI{"debugging bypass gsi", false, this};
};

// -------------------------------------------------------------------------------------------------

class Interpolation : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::gsi::Interpolation";}

  typedef InterpolationParameters Parameters_;

  Interpolation(const oops::GeometryData &,
                const std::vector<size_t> &,
                const oops::Variables &,
                const Parameters_ &,
                const atlas::FieldSet &,
                const atlas::FieldSet &,
                const std::vector<atlas::FieldSet> &);
  virtual ~Interpolation();

  const oops::GeometryData & innerGeometryData() const override {return *innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<const oops::GeometryData> innerGeometryData_;
  oops::Variables innerVars_;

  // Interpolation object
  std::unique_ptr<UnstructuredInterpolation> interpolator_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
