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

#include <boost/noncopyable.hpp>

#include "oops/base/GeometryData.h"
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

class SaberOuterBlockBase : public util::Printable, private boost::noncopyable {
 public:
  SaberOuterBlockBase() {}
  virtual ~SaberOuterBlockBase() {}

  virtual const oops::GeometryData & inputGeometryData() const = 0;
  virtual const oops::Variables & inputVars() const = 0;

  virtual void multiply(atlas::FieldSet &) const = 0;
  virtual void multiplyAD(atlas::FieldSet &) const = 0;
  virtual void calibrationInverseMultiply(atlas::FieldSet &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

class SaberOuterBlockFactory;

// -----------------------------------------------------------------------------

class SaberOuterBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberOuterBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberOuterBlockFactory>
    saberOuterBlockParameters{"saber block name", this};
};

// =============================================================================

class SaberOuterBlockFactory {
 public:
  static SaberOuterBlockBase * create(const oops::GeometryData &,
                                      const std::vector<size_t> &,
                                      const oops::Variables &,
                                      const SaberBlockParametersBase &,
                                      const atlas::FieldSet &,
                                      const atlas::FieldSet &,
                                      const std::vector<atlas::FieldSet> &);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberOuterBlockFactory() = default;

 protected:
  explicit SaberOuterBlockFactory(const std::string &name);

 private:
  virtual SaberOuterBlockBase * make(const oops::GeometryData &,
                                     const std::vector<size_t> &,
                                     const oops::Variables &,
                                     const SaberBlockParametersBase &,
                                     const atlas::FieldSet &,
                                     const atlas::FieldSet &,
                                     const std::vector<atlas::FieldSet> &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberOuterBlockFactory * > & getMakers() {
    static std::map < std::string, SaberOuterBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberOuterBlockMaker : public SaberOuterBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  SaberOuterBlockBase * make(const oops::GeometryData & outputGeometryData,
                             const std::vector<size_t> & activeVariableSizes,
                             const oops::Variables & outputVars,
                             const SaberBlockParametersBase & params,
                             const atlas::FieldSet & xb,
                             const atlas::FieldSet & fg,
                             const std::vector<atlas::FieldSet> & fsetVec) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(outputGeometryData, activeVariableSizes,
                 outputVars, stronglyTypedParams, xb, fg, fsetVec);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberOuterBlockMaker(const std::string & name) : SaberOuterBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace saber
