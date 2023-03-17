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
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

class SaberCentralBlockBase : public util::Printable, private boost::noncopyable {
 public:
  SaberCentralBlockBase() {}
  virtual ~SaberCentralBlockBase() {}

  virtual void randomize(atlas::FieldSet &) const = 0;
  virtual void multiply(atlas::FieldSet &) const = 0;

  void adjointTest(const eckit::mpi::Comm &,
                   const oops::GeometryData &,
                   const std::vector<size_t> &,
                   const oops::Variables &,
                   const double & adjointTolerance = 1.0e-12) const;

 private:
  virtual void print(std::ostream &) const = 0;
};

// =============================================================================

class SaberCentralBlockFactory;

// -----------------------------------------------------------------------------

class SaberCentralBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberCentralBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberCentralBlockFactory>
    saberCentralBlockParameters{"saber block name", this};
};

// =============================================================================

class SaberCentralBlockFactory {
 public:
  static SaberCentralBlockBase * create(const oops::GeometryData &,
                                        const std::vector<size_t> &,
                                        const oops::Variables &,
                                        const SaberBlockParametersBase &,
                                        const atlas::FieldSet &,
                                        const atlas::FieldSet &,
                                        const std::vector<atlas::FieldSet> &,
                                        const size_t & timeRank = 0);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberCentralBlockFactory() = default;

 protected:
  explicit SaberCentralBlockFactory(const std::string &name);

 private:
  virtual SaberCentralBlockBase * make(const oops::GeometryData &,
                                       const std::vector<size_t> &,
                                       const oops::Variables &,
                                       const SaberBlockParametersBase &,
                                       const atlas::FieldSet &,
                                       const atlas::FieldSet &,
                                       const std::vector<atlas::FieldSet> &,
                                       const size_t &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberCentralBlockFactory * > & getMakers() {
    static std::map < std::string, SaberCentralBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberCentralBlockMaker : public SaberCentralBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  SaberCentralBlockBase * make(const oops::GeometryData & geometryData,
                               const std::vector<size_t> & activeVariableSizes,
                               const oops::Variables & outerVars,
                               const SaberBlockParametersBase & params,
                               const atlas::FieldSet & xb,
                               const atlas::FieldSet & fg,
                               const std::vector<atlas::FieldSet> & fsetVec,
                               const size_t & timeRank) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(geometryData, activeVariableSizes,
                 outerVars, stronglyTypedParams, xb, fg, fsetVec, timeRank);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberCentralBlockMaker(const std::string & name) : SaberCentralBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace saber
