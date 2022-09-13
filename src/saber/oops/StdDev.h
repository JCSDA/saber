/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_STDDEV_H_
#define SABER_OOPS_STDDEV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class StdDevParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(StdDevParameters, SaberBlockParametersBase)
};

// -----------------------------------------------------------------------------

class StdDev : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::StdDev";}

  typedef StdDevParameters Parameters_;

  StdDev(const eckit::mpi::Comm &,
         const atlas::FunctionSpace &,
         const atlas::FieldSet &,
         const std::vector<size_t> &,
         const Parameters_ &,
         const atlas::FieldSet &,
         const atlas::FieldSet &,
         const std::vector<atlas::FieldSet> &);
  virtual ~StdDev();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  atlas::FieldSet stdDevFset_;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_STDDEV_H_
