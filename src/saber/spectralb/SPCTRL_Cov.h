/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"

#include "saber/spectralb/spectralb.h"
#include "saber/spectralb/spectralbParameters.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

class SPCTRL_COVParameters : public SaberCentralBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SPCTRL_COVParameters, SaberCentralBlockParametersBase)

 public:
  oops::RequiredParameter<spectralbParameters> spectralbParams{"spectralb", this};
};

// -----------------------------------------------------------------------------

class SPCTRL_COV : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::lfricspectralb::SPCTRL_COV";}

  typedef SPCTRL_COVParameters Parameters_;

  SPCTRL_COV(const eckit::mpi::Comm &,
     const atlas::FunctionSpace &,
     const atlas::FieldSet &,
     const std::vector<size_t> &,
     const eckit::Configuration &,
     const atlas::FieldSet &,
     const atlas::FieldSet &,
     const std::vector<atlas::FieldSet> &);
  virtual ~SPCTRL_COV();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SpectralB> spectralb_;
};

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
