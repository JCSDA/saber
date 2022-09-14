/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_BUMP_BUMP_STDDEV_H_
#define SABER_BUMP_BUMP_STDDEV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"

#include "saber/bump/BUMP.h"
#include "saber/bump/BUMP_Parameters.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class BUMP_StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BUMP_StdDevParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<BUMP_Parameters> bumpParams{"bump", this};
};

// -----------------------------------------------------------------------------


class BUMP_StdDev : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::BUMP_StdDev";}

  typedef BUMP_StdDevParameters Parameters_;

  BUMP_StdDev(const eckit::mpi::Comm &,
              const atlas::FunctionSpace &,
              const atlas::FieldSet &,
              const std::vector<size_t> &,
              const Parameters_ &,
              const atlas::FieldSet &,
              const atlas::FieldSet &,
              const std::vector<atlas::FieldSet> &);
  virtual ~BUMP_StdDev();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<BUMP> bump_;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_BUMP_BUMP_STDDEV_H_
