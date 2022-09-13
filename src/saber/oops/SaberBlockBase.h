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

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include <boost/noncopyable.hpp>

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

class SaberBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberBlockBase(const SaberBlockParametersBase & params);
  virtual ~SaberBlockBase() {}

  virtual void randomize(atlas::FieldSet &) const = 0;
  virtual void multiply(atlas::FieldSet &) const = 0;
  virtual void inverseMultiply(atlas::FieldSet &) const = 0;
  virtual void multiplyAD(atlas::FieldSet &) const = 0;
  virtual void inverseMultiplyAD(atlas::FieldSet &) const = 0;

  bool iterativeInverse() const {return iterativeInverse_;}
  const std::string name() const {return name_;}
 private:
  virtual void print(std::ostream &) const = 0;
  bool iterativeInverse_;
  std::string name_;
};

// =============================================================================

class SaberBlockFactory;

// -----------------------------------------------------------------------------

class SaberBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberBlockFactory>
    saberBlockParameters{"saber block name", this};
};

// =============================================================================

class SaberBlockFactory {

 public:
  static SaberBlockBase * create(const atlas::FunctionSpace &,
                                 const atlas::FieldSet &,
                                 const std::vector<size_t> &,
                                 const SaberBlockParametersBase &,
                                 const atlas::FieldSet &,
                                 const atlas::FieldSet &,
                                 const std::vector<atlas::FieldSet> &);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberBlockFactory() = default;

 protected:
  explicit SaberBlockFactory(const std::string &name);

 private:
  virtual SaberBlockBase * make(const atlas::FunctionSpace &,
                                const atlas::FieldSet &,
                                const std::vector<size_t> &,
                                const SaberBlockParametersBase &,
                                const atlas::FieldSet &,
                                const atlas::FieldSet &,
                                const std::vector<atlas::FieldSet> &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberBlockFactory * > & getMakers() {
    static std::map < std::string, SaberBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberBlockMaker : public SaberBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  SaberBlockBase * make(const atlas::FunctionSpace & functionSpace,
                        const atlas::FieldSet & extraFields,
                        const std::vector<size_t> & variableSizes,
                        const SaberBlockParametersBase & params,
                        const atlas::FieldSet & xb,
                        const atlas::FieldSet & fg,
                        const std::vector<atlas::FieldSet> & fsetVec) override {
    return new T(functionSpace, extraFields, variableSizes,
                 dynamic_cast<const Parameters_&>(params), xb, fg, fsetVec);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberBlockMaker(const std::string & name) : SaberBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_SABERBLOCKBASE_H_
