/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"

#include "saber/bifourier/BifourierTransformBase.h"
#include "saber/bifourier/BifourierTransformStore.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierIDParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BifourierIDParameters, SaberBlockParametersBase)

  oops::Variables mandatoryActiveVars() const override
    {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class BifourierID : public SaberCentralBlockBase {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierID";}

  typedef BifourierIDParameters Parameters_;

  BifourierID(const oops::GeometryData &,
              const oops::Variables &,
              const eckit::Configuration &,
              const BifourierIDParameters &,
              const oops::FieldSet3D &,
              const oops::FieldSet3D &);
  virtual ~BifourierID();

  size_t ctlVecSize() const override
    {return trans_->ctlVecSize();}
  void randomCtlVec(atlas::Field & cv,
                    const size_t & offset) const override
    {trans_->randomCtlVec(cv, centralVars(), offset);}
  void multiplySqrt(const atlas::Field &,
                    oops::FieldSet3D &,
                    const size_t &) const override;
  void multiplySqrtAD(const oops::FieldSet3D &,
                      atlas::Field &,
                      const size_t &) const override;

  void read() override
    {}

 private:
  // Communicator
  const eckit::mpi::Comm & comm_;

  // Spectral transform
  const BifourierTransformStore transStore_;
  const std::shared_ptr<BifourierTransformBase> trans_;

  // Private methods

  // Print
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
