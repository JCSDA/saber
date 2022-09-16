/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

#include "saber/oops/SaberOuterBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberOuterBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberOuterBlockBase(const eckit::Configuration & conf);
  virtual ~SaberOuterBlockBase() {}

  const oops::Variables inputVars() {return inputVars_;}
  const atlas::FunctionSpace inputFunctionSpace() {return inputFunctionSpace_;}
  const atlas::FieldSet inputExtraFields() {return inputExtraFields_;}

  virtual void multiply(atlas::FieldSet &) const = 0;
  virtual void multiplyAD(atlas::FieldSet &) const = 0;
  virtual void calibrationInverseMultiply(atlas::FieldSet &) const = 0;

  const std::string name() const {return name_;}

 protected:
  oops::Variables inputVars_;
  atlas::FunctionSpace inputFunctionSpace_;
  atlas::FieldSet inputExtraFields_;

 private:
  virtual void print(std::ostream &) const = 0;
  std::string name_;
};

// =============================================================================

class SaberOuterBlockFactory;

// -----------------------------------------------------------------------------

class SaberOuterBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberOuterBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberOuterBlockParametersBase, SaberOuterBlockFactory>
    saberOuterBlockParameters{"saber block name", this};
};

// =============================================================================

class SaberOuterBlockFactory {

 public:
  static SaberOuterBlockBase * create(const eckit::mpi::Comm &,
                                      const atlas::FunctionSpace &,
                                      const atlas::FieldSet &,
                                      const std::vector<size_t> &,
                                      const eckit::Configuration &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &,
                                      const std::vector<atlas::FieldSet> &);

  static std::unique_ptr<SaberOuterBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberOuterBlockFactory() = default;

 protected:
  explicit SaberOuterBlockFactory(const std::string &name);

 private:
  virtual SaberOuterBlockBase * make(const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const atlas::FieldSet &,
                                     const std::vector<size_t> &,
                                     const eckit::Configuration &,
                                     const atlas::FieldSet &,
                                     const atlas::FieldSet &,
                                     const std::vector<atlas::FieldSet> &) = 0;

  virtual std::unique_ptr<SaberOuterBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberOuterBlockFactory * > & getMakers() {
    static std::map < std::string, SaberOuterBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberOuterBlockMaker : public SaberOuterBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  SaberOuterBlockBase * make(const eckit::mpi::Comm & comm,
                             const atlas::FunctionSpace & outputFunctionSpace,
                             const atlas::FieldSet & outputExtraFields,
                             const std::vector<size_t> & activeVariableSizes,
                             const eckit::Configuration & conf,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const std::vector<atlas::FieldSet> & fsetVec) override {
    return new T(comm, outputFunctionSpace, outputExtraFields, activeVariableSizes,
                 conf, xb, fg, fsetVec);
  }

  std::unique_ptr<SaberOuterBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberOuterBlockMaker(const std::string & name) : SaberOuterBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace saber
