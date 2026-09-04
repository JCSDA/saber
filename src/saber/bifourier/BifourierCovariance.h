/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>

#include "atlas/field.h"

#include "oops/base/FieldSets.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variable.h"

#include "saber/bifourier/BifourierCovarianceImpl.h"
#include "saber/blocks/SaberCentralBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierCovariance : public SaberCentralBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierCovariance";}

  typedef BifourierCovarianceImplParameters Parameters_;

  BifourierCovariance(const oops::GeometryData &,
                      const oops::Variables &,
                      const eckit::Configuration &,
                      const Parameters_ &,
                      const oops::FieldSet3D &,
                      const oops::FieldSet3D &);
  virtual ~BifourierCovariance();

  size_t ctlVecSize() const override
    {return covar_->ctlVecSize();}
  void randomCtlVec(atlas::Field & cv,
                    const size_t & offset) const override
    {covar_->randomCtlVec(cv, offset);}
  void multiplySqrt(const atlas::Field & cv,
                    oops::FieldSet3D & fset,
                    const size_t & offset) const override
    {covar_->multiplySqrt(cv, fset, offset);}
  void multiplySqrtAD(const oops::FieldSet3D & fset,
                      atlas::Field & cv,
                      const size_t & offset) const override
    {covar_->multiplySqrtAD(fset, cv, offset);}

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
