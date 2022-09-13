/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_SPECTRALB_SPCTRL_COV_H_
#define SABER_SPECTRALB_SPCTRL_COV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

#include "saber/spectralb/spectralb.h"
#include "saber/spectralb/spectralbParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------
class SPCTRL_COVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SPCTRL_COVParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<spectralbParameters> spectralbParams{"spectralb", this};
};

// -----------------------------------------------------------------------------

class SPCTRL_COV : public SaberBlockBase {

 public:
  static const std::string classname() {return "saber::lfricspectralb::SPCTRL_COV";}

  typedef SPCTRL_COVParameters Parameters_;

  SPCTRL_COV(const eckit::mpi::Comm &,
             const atlas::FunctionSpace &,
             const atlas::FieldSet &,
             const std::vector<size_t> &,
             const Parameters_ &,
             const atlas::FieldSet &,
             const atlas::FieldSet &,
             const std::vector<atlas::FieldSet> &);
  virtual ~SPCTRL_COV();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SpectralB> spectralb_;
};

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_SPECTRALB_SPCTRL_COV_H_
