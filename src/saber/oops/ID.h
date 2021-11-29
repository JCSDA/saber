/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ID_H_
#define SABER_OOPS_ID_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

class IDParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(IDParameters, SaberBlockParametersBase)
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class ID : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL>  Geometry_;

 public:
  static const std::string classname() {return "saber::ID";}

  typedef IDParameters Parameters_;

  ID(const Geometry_ &, const Parameters_ &);
  virtual ~ID();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ID<MODEL>::ID(const Geometry_ & resol, const Parameters_ & params)
  : SaberBlockBase<MODEL>(params)
{
  oops::Log::trace() << classname() << "::ID starting" << std::endl;
  oops::Log::trace() << classname() << "::ID done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ID<MODEL>::~ID() {
  oops::Log::trace() << classname() << "::~ID starting" << std::endl;
  util::Timer timer(classname(), "~ID");
  oops::Log::trace() << classname() << "::~ID done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::randomize(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::multiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::inverseMultiply(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::multiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::inverseMultiplyAD(atlas::FieldSet * atlasFieldSet) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ID<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ID_H_
