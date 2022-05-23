/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SPCTRL_COV_H_
#define SABER_OOPS_SPCTRL_COV_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/spectralb/spectralb.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------
template <typename MODEL>
class SPCTRL_COVParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SPCTRL_COVParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<spectralbParameters<MODEL>> spectralbParams{"spectralb", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class SPCTRL_COV : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;
  typedef SpectralB<MODEL>      SpectralB_;

 public:
  static const std::string classname() {return "saber::lfricspectralb::SPCTRL_COV";}

  typedef SPCTRL_COVParameters<MODEL> Parameters_;

  SPCTRL_COV(const Geometry_ &,
                  const Parameters_ &,
                  const State_ &,
                  const State_ &);
  virtual ~SPCTRL_COV();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SpectralB_> spectralb_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
SPCTRL_COV<MODEL>::SPCTRL_COV(const Geometry_ & resol, const Parameters_ & params,
                                        const State_ & xb, const State_ & fg)
  : SaberBlockBase<MODEL>(params), spectralb_()
{
  oops::Log::trace() << classname() << "::SPCTRL_COV starting" << std::endl;

  // Setup and check input/ouput variables
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  ASSERT(inputVars == outputVars);

  // Active variables
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= inputVars);
  } else {
    activeVars += inputVars;
  }

  // Initialize SpectralB_
  spectralb_.reset(new SpectralB_(resol, activeVars, params.spectralbParams.value()));

  oops::Log::trace() << classname() << "::SPCTRL_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
SPCTRL_COV<MODEL>::~SPCTRL_COV() {
  oops::Log::trace() << classname() << "::~SPCTRL_COV starting" << std::endl;
  util::Timer timer(classname(), "~SPCTRL_COV");
  oops::Log::trace() << classname() << "::~SPCTRL_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  ABORT("SPCTRL_COV<MODEL>::randomize: not implemented");
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  spectralb_->multiply_InterpAndCov(atlasFieldSet);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  ABORT("SPCTRL_COV<MODEL>::inverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  ABORT("SPCTRL_COV<MODEL>::multiplyAD: not implemented");
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("SPCTRL_COV<MODEL>::inverseMultiplyAD: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void SPCTRL_COV<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber

#endif  // SABER_OOPS_SPCTRL_COV_H_
