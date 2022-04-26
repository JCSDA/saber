/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_SABERBLOCKBASE_H_
#define SABER_OOPS_SABERBLOCKBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberBlockBase(const SaberBlockParametersBase & params);
  virtual ~SaberBlockBase() {}

  virtual void randomize(atlas::FieldSet *) const = 0;
  virtual void multiply(atlas::FieldSet *) const = 0;
  virtual void inverseMultiply(atlas::FieldSet *) const = 0;
  virtual void multiplyAD(atlas::FieldSet *) const = 0;
  virtual void inverseMultiplyAD(atlas::FieldSet *) const = 0;

  bool iterativeInverse() const {return iterativeInverse_;}
 private:
  virtual void print(std::ostream &) const = 0;
  bool iterativeInverse_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SaberBlockBase<MODEL>::SaberBlockBase(const SaberBlockParametersBase & params)
  : iterativeInverse_(params.iterativeInverse.value()) {}

// =============================================================================

template <typename MODEL>
class SaberBlockFactory;

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberBlockFactory<MODEL>>
    saberBlockParameters{"saber block name", this};
};

// =============================================================================

template <typename MODEL>
class SaberBlockFactory {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;

 public:
  static SaberBlockBase<MODEL> * create(const Geometry_ &,
                                        const SaberBlockParametersBase &,
                                        const State_ & xb,
                                        const State_ & fg);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberBlockFactory() = default;

 protected:
  explicit SaberBlockFactory(const std::string &name);

 private:
  virtual SaberBlockBase<MODEL> * make(const Geometry_ &,
                                       const SaberBlockParametersBase &,
                                       const State_ &,
                                       const State_ &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberBlockFactory<MODEL> * > & getMakers() {
    static std::map < std::string, SaberBlockFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class SaberBlockMaker : public SaberBlockFactory<MODEL> {
  typedef typename T::Parameters_ Parameters_;
  typedef oops::Geometry<MODEL>   Geometry_;
  typedef oops::State<MODEL>      State_;

  SaberBlockBase<MODEL> * make(const Geometry_ & geom,
                               const SaberBlockParametersBase & params,
                               const State_ & xb,
                               const State_ & fg) override {
    return new T(geom, dynamic_cast<const Parameters_&>(params), xb, fg);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberBlockMaker(const std::string & name) : SaberBlockFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SaberBlockFactory<MODEL>::SaberBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberBlockFactory." << std::endl;
    ABORT("Element already registered in saber::SaberBlockFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
SaberBlockBase<MODEL> * SaberBlockFactory<MODEL>::create(const Geometry_ & geom,
                                                         const SaberBlockParametersBase & params,
                                                         const State_& xb,
                                                         const State_ & fg) {
  oops::Log::trace() << "SaberBlockBase<MODEL>::create starting" << std::endl;
  const std::string &id = params.saberBlockName.value().value();
  typename std::map<std::string, SaberBlockFactory<MODEL>*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberBlockFactory." << std::endl;
    ABORT("Element does not exist in saber::SaberBlockFactory.");
  }
  SaberBlockBase<MODEL> * ptr = jsb->second->make(geom, params, xb, fg);
  oops::Log::trace() << "SaberBlockBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<SaberBlockParametersBase>
SaberBlockFactory<MODEL>::createParameters(const std::string &name) {
  typename std::map<std::string, SaberBlockFactory<MODEL>*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in saber::SaberBlockFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKBASE_H_
