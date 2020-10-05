/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCEID_H_
#define SABER_OOPS_ERRORCOVARIANCEID_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -------------------------------------------------------------------------------------------------

/// Identity model space error covariance

template <typename MODEL>
class ErrorCovarianceID : public oops::ModelSpaceCovarianceBase<MODEL>,
                          public util::Printable,
                          private util::ObjectCounter<ErrorCovarianceID<MODEL>> {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment<MODEL>   Increment_;
  typedef oops::State<MODEL>       State_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceID";}

  ErrorCovarianceID(const Geometry_ &, const oops::Variables &, const eckit::Configuration &,
                    const State_ &, const State_ &);
  virtual ~ErrorCovarianceID();

 private:
  ErrorCovarianceID(const ErrorCovarianceID&);
  ErrorCovarianceID& operator=(const ErrorCovarianceID&);

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;
  void print(std::ostream &) const override;
};

// =================================================================================================

template<typename MODEL>
ErrorCovarianceID<MODEL>::ErrorCovarianceID(const Geometry_ & resol, const oops::Variables & vars,
                                            const eckit::Configuration & conf, const State_ & xb,
                                            const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf)
{
  oops::Log::trace() << "ErrorCovarianceID::ErrorCovarianceID starting" << std::endl;
  oops::Log::trace() << "ErrorCovarianceID::ErrorCovarianceID done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
ErrorCovarianceID<MODEL>::~ErrorCovarianceID() {
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::~ErrorCovarianceID starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovarianceID");
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::~ErrorCovarianceID done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceID<MODEL>::doRandomize(Increment_ & dx) const {
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  dx.random();
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doRandomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceID<MODEL>::doMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  dxo = dxi;
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doMultiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceID<MODEL>::doInverseMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  dxo = dxi;
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::doInverseMultiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceID<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovarianceID<MODEL> B = I";
  oops::Log::trace() << "ErrorCovarianceID<MODEL>::print done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCEID_H_
