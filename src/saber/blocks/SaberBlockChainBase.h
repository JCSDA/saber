/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/Logger.h"

namespace atlas {
  class Field;
  class FieldSet;
  class FunctionSpace;
}

namespace oops {
  class FieldSet3D;
  class FieldSet4D;
  class FieldSets;
  template <class MODEL> class Geometry;
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// Base class for SABER block chains that have a self-adjoint central block
/// (ensemble and non-ensemble).
class SaberBlockChainBase {
 public:
  SaberBlockChainBase() = default;
  virtual ~SaberBlockChainBase() = default;

  virtual void randomize(oops::FieldSet4D &) const = 0;
  virtual void multiply(oops::FieldSet4D &) const = 0;
  virtual size_t ctlVecSize() const = 0;
  virtual void randomCtlVec(atlas::Field &, const size_t &) const = 0;
  virtual void multiplySqrt(const atlas::Field &, oops::FieldSet4D &, const size_t &) const = 0;
  virtual void multiplySqrtAD(const oops::FieldSet4D &, atlas::Field &, const size_t &)
    const = 0;
  virtual const atlas::FunctionSpace & outerFunctionSpace() const = 0;
  virtual const oops::Variables & outerVariables() const = 0;

  /// @brief Diagonal variance of this block chain on its outer function space.
  virtual oops::FieldSet3D variance() const = 0;
};

template<typename MODEL>
class SaberBlockChainFactory {
 public:
  typedef oops::Geometry<MODEL> Geometry_;

  static std::unique_ptr<SaberBlockChainBase> create(const Geometry_ &,
                                                     const oops::Variables &,
                                                     oops::FieldSet4D &,
                                                     oops::FieldSet4D &,
                                                     const eckit::Configuration &);

  virtual ~SaberBlockChainFactory() = default;

 protected:
  explicit SaberBlockChainFactory(const std::string &);

 private:
  virtual std::unique_ptr<SaberBlockChainBase> make(const Geometry_ &,
                                                    const oops::Variables &,
                                                    oops::FieldSet4D &,
                                                    oops::FieldSet4D &,
                                                    const eckit::Configuration &) = 0;

  static std::map <std::string, SaberBlockChainFactory<MODEL> *> & getMakers() {
    static std::map <std::string, SaberBlockChainFactory<MODEL> *> makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class SaberBlockChainMaker : public SaberBlockChainFactory<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;

  std::unique_ptr<SaberBlockChainBase> make(const Geometry_ & geom,
                                            const oops::Variables & outerVars,
                                            oops::FieldSet4D & fset4dXb,
                                            oops::FieldSet4D & fset4dFg,
                                            const eckit::Configuration & conf) override {
    return std::make_unique<T>(geom, outerVars, fset4dXb, fset4dFg, conf);
  }

 public:
  explicit SaberBlockChainMaker(const std::string & name) : SaberBlockChainFactory<MODEL>(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
SaberBlockChainFactory<MODEL>::SaberBlockChainFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name + " already registered in saber::SaberBlockChainFactory.",
                              Here());
  getMakers()[name] = this;
}

template <typename MODEL>
std::unique_ptr<SaberBlockChainBase>
SaberBlockChainFactory<MODEL>::create(const Geometry_ & geom,
                                      const oops::Variables & outerVars,
                                      oops::FieldSet4D & fset4dXb,
                                      oops::FieldSet4D & fset4dFg,
                                      const eckit::Configuration & conf) {
  oops::Log::trace() << "SaberBlockChainFactory<MODEL>::create starting" << std::endl;
  std::string name = "parametric";
  if (conf.has("covariance type")) {
    name = conf.getString("covariance type");
  }
  typename std::map<std::string, SaberBlockChainFactory<MODEL>*>::iterator jbc =
    getMakers().find(name);
  if (jbc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in saber::SaberBlockChainFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<SaberBlockChainBase> ptr =
    jbc->second->make(geom, outerVars, fset4dXb, fset4dFg, conf);
  oops::Log::trace() << "SaberBlockChainFactory<MODEL>::create done" << std::endl;
  return ptr;
}

}  // namespace saber
