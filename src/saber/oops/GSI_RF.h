/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_GSI_RF_H_
#define SABER_OOPS_GSI_RF_H_

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"

#include "saber/oops/GSI.h"
#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class GSI_RFParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(GSI_RFParameters, SaberBlockParametersBase)

 public:
  oops::RequiredParameter<GSI_Parameters> gsiParams{"gsi", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class GSI_RF : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef GSI<MODEL>            GSI_;

 public:
  static const std::string classname() {return "saber::GSI_RF";}

  typedef GSI_RFParameters Parameters_;

  GSI_RF(const Geometry_ &, const Parameters_ &);
  virtual ~GSI_RF();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<GSI_> gsi_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
GSI_RF<MODEL>::GSI_RF(const Geometry_ & resol, const Parameters_ & params)
  : SaberBlockBase<MODEL>(params)
{
  oops::Log::trace() << classname() << "::GSI_RF starting" << std::endl;

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

  // Initialize GSI
  gsi_.reset(new GSI_(resol, activeVars, params.gsiParams.value()));

  oops::Log::trace() << classname() << "::GSI_RF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
GSI_RF<MODEL>::~GSI_RF() {
  oops::Log::trace() << classname() << "::~GSI_RF starting" << std::endl;
  util::Timer timer(classname(), "~GSI_RF");
  oops::Log::trace() << classname() << "::~GSI_RF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void GSI_RF<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_GSI_RF_H_
