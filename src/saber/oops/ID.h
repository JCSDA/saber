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

class ID : public SaberBlockBase {
 public:
  static const std::string classname() {return "saber::ID";}

  typedef IDParameters Parameters_;

  ID(const atlas::FunctionSpace &,
     const atlas::FieldSet &,
     const Parameters_ &,
     const atlas::FieldSet &,
     const atlas::FieldSet &);
  virtual ~ID();

  void initialize(const std::vector<atlas::FieldSet> &) override;

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;
  void inverseMultiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void inverseMultiplyAD(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

ID::ID(const atlas::FunctionSpace & functionSpace,
       const atlas::FieldSet & extraFields,
       const Parameters_ & params,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg)
  : SaberBlockBase(params)
{
  oops::Log::trace() << classname() << "::ID starting" << std::endl;
  oops::Log::trace() << classname() << "::ID done" << std::endl;
}

// -----------------------------------------------------------------------------

ID::~ID() {
  oops::Log::trace() << classname() << "::~ID starting" << std::endl;
  util::Timer timer(classname(), "~ID");
  oops::Log::trace() << classname() << "::~ID done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::initialize(const std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << classname() << "::initialize starting" << std::endl;
  oops::Log::trace() << classname() << "::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ID::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ID_H_
