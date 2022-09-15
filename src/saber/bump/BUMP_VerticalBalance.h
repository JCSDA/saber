/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_BUMP_VERTICALBALANCE_H_
#define SABER_BUMP_BUMP_VERTICALBALANCE_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class BUMP_VerticalBalanceParameters : public SaberOuterBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_VerticalBalanceParameters, SaberOuterBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------

class BUMP_VerticalBalance : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::BUMP_VerticalBalance";}

  typedef BUMP_VerticalBalanceParameters Parameters_;

  BUMP_VerticalBalance(const eckit::mpi::Comm &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const eckit::Configuration &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);
  virtual ~BUMP_VerticalBalance();

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void calibrationInverseMultiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_BUMP_BUMP_VERTICALBALANCE_H_
