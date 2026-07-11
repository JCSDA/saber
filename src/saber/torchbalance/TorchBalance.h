/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
* (C) Copyright 2025- UCAR
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <map>
#include <string>
#include <vector>

#include "oops/base/Variables.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"


// --------------------------------------------------------------------------------------

namespace saber {

// --------------------------------------------------------------------------------------

class SurfaceEmulatorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SurfaceEmulatorParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<std::string> path{"path", this};
  oops::RequiredParameter<std::vector<std::string>> jacobianWrt{"jacobian wrt", this};
  oops::OptionalParameter<std::string> maskVariable{"jacobian masking.variable", this};
  oops::OptionalParameter<int> maskLevel{"jacobian masking.level", this};
};

class VerticalEmulatorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VerticalEmulatorParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<std::string> path{"path", this};
  oops::RequiredParameter<std::vector<std::string>> jacobianWrt{"jacobian wrt", this};
  oops::OptionalParameter<std::string> maskVariable{"jacobian masking.variable", this};
  oops::OptionalParameter<int> maskLevel{"jacobian masking.level", this};
};

class TorchBalanceParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(TorchBalanceParameters, saber::SaberBlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {
    // No mandatory variables - fully generic and configuration-driven
    return oops::Variables();
  }

  oops::Parameter<oops::Variables> innerVars{"inner variables", oops::Variables(), this};
  oops::RequiredParameter<std::vector<SurfaceEmulatorParameters>> surfaceEmulators{
    "surface emulators", this};
  oops::OptionalParameter<std::vector<VerticalEmulatorParameters>> verticalEmulators{
    "vertical emulators", this};
  oops::Parameter<bool> saveJacobians{"save jacobians", false, this};
};

// --------------------------------------------------------------------------------------

class TorchBalance : public SaberOuterBlockBase {
 public:
  static const std::string classname() { return "saber::TorchBalance"; }
  typedef TorchBalanceParameters Parameters_;

  TorchBalance(const oops::GeometryData &,
               const oops::Variables &,
               const eckit::Configuration &,
               const Parameters_ &,
               const oops::FieldSet3D &,
               const oops::FieldSet3D &);
  virtual ~TorchBalance() = default;

  const oops::Variables & innerVars() const override {return innerVars_;}
  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

  /// Parse a Jacobian field name of the form "d{numerator}_div_d{denominator}"
  /// and return {numerator, denominator}. Returns false if the name is not valid.
  static bool parseJacobianName(const std::string & jacName,
                                std::string & numerator,
                                std::string & denominator);

  /// Look up the level metadata for a given Jacobian field name.
  /// Sets all four level fields and returns true on success.
  /// Returns false if the Jacobian name is not found in jacLevelMetadata_.
  bool getJacobianLevels(const std::string & jacName,
                         int & denominator_level,
                         int & numerator_level,
                         int & nInputLevels,
                         int & nOutputLevels) const;

 private:
  void print(std::ostream &) const override;
  oops::Variables innerVars_;
  const oops::GeometryData & innerGeometryData_;
  atlas::FieldSet jac_;

  // Store level metadata for each Jacobian field.
  // Surface emulators: nInputLevels = nOutputLevels = 1, specific level indices stored.
  // Vertical emulators: denominator_level = numerator_level = 0, compact vertical block used.
  struct JacobianLevels {
    int denominator_level;
    int numerator_level;
    int nInputLevels  = 1;
    int nOutputLevels = 1;
  };
  std::map<std::string, JacobianLevels> jacLevelMetadata_;
};

// --------------------------------------------------------------------------------------

}  // namespace saber
