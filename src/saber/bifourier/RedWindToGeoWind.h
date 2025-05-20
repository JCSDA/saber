/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameters.h"

#include "saber/bifourier/BiperiodizationImpl.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class RedWindToGeoWindParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(RedWindToGeoWindParameters, SaberBlockParametersBase)

 public:
  // Biperiodization parameters
  oops::OptionalParameter<BiperiodizationImplParameters> biperParams{"biperiodization",
    this};

  // Outer spherical winds
  oops::Parameter<bool> outerSphericalWinds{"outer spherical winds", false, this};

  // Output file
  oops::OptionalParameter<eckit::LocalConfiguration> outputFile{"output file", this};

  oops::Variables mandatoryActiveVars() const override {
    oops::Variables mandatoryActiveVars;
    mandatoryActiveVars.push_back("reduced_x_wind");
    mandatoryActiveVars.push_back("reduced_y_wind");
    if (outerSphericalWinds) {
      mandatoryActiveVars.push_back("eastward_wind");
      mandatoryActiveVars.push_back("northward_wind");
    } else {
      mandatoryActiveVars.push_back("geographical_x_wind");
      mandatoryActiveVars.push_back("geographical_y_wind");
    }
    return mandatoryActiveVars;
  }
};

// -----------------------------------------------------------------------------

class RedWindToGeoWind : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::RedWindToGeoWind";}

  typedef RedWindToGeoWindParameters Parameters_;

  RedWindToGeoWind(const oops::GeometryData &,
                   const oops::Variables &,
                   const eckit::Configuration &,
                   const Parameters_ &,
                   const oops::FieldSet3D &,
                   const oops::FieldSet3D &);
  virtual ~RedWindToGeoWind() = default;

  const oops::GeometryData & innerGeometryData() const override
    {return innerGeometryData_;}
  const oops::Variables & innerVars() const override
    {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const;

 private:
  // Inner geometry data
  const oops::GeometryData & innerGeometryData_;

  // Communicator
  const eckit::mpi::Comm & comm_;

  // Inner variables
  oops::Variables innerVars_;

  // Parameters
  Parameters_ params_;

  // Map factor and Jacobian coefficients FieldSet
  oops::FieldSet3D data_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
