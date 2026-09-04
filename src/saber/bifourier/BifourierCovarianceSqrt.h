/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variable.h"

#include "saber/bifourier/BifourierCovarianceImpl.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierCovarianceSqrt : public SaberOuterBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierCovarianceSqrt";}

  typedef BifourierCovarianceImplParameters Parameters_;

  BifourierCovarianceSqrt(const oops::GeometryData &,
                          const oops::Variables &,
                          const eckit::Configuration &,
                          const Parameters_ &,
                          const oops::FieldSet3D &,
                          const oops::FieldSet3D &);
  virtual ~BifourierCovarianceSqrt();

  const oops::GeometryData & innerGeometryData() const override
    {return covar_->innerGeometryData();}
  const oops::Variables & innerVars() const override
    {return covar_->innerVars();}

  void multiply(oops::FieldSet3D & fset) const override
    {covar_->multiplySqrt(fset);}
  void multiplyAD(oops::FieldSet3D & fset) const override
    {covar_->multiplySqrtAD(fset);}
  void leftInverseMultiply(oops::FieldSet3D & fset) const override
    {covar_->leftInverseMultiply(fset);}
  void read() override
    {covar_->read();}

  void directCalibration(const oops::FieldSets & fsetEns) override
    {covar_->directCalibration(fsetEns);}

  void iterativeCalibrationInit() override
    {covar_->iterativeCalibrationInit();}
  void iterativeCalibrationUpdate(const oops::FieldSet3D & fset) override
    {covar_->iterativeCalibrationUpdate(fset);}
  void iterativeCalibrationFinal() override
    {covar_->iterativeCalibrationFinal();}

  void write() const override
    {covar_->write();}

 protected:
  // Covariance implementation
  std::unique_ptr<BifourierCovarianceImpl> covar_;

  // Print
  void print(std::ostream & os) const override
    {covar_->print(os);}
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
