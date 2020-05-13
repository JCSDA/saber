/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCEBUMP_H_
#define SABER_OOPS_ERRORCOVARIANCEBUMP_H_

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

#include "saber/oops/OoBump.h"
#include "saber/oops/ParametersBUMP.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

/// Model space error covariance with BUMP

template <typename MODEL>
class ErrorCovarianceBUMP : public oops::ModelSpaceCovarianceBase<MODEL>,
                            public util::Printable,
                            private util::ObjectCounter<ErrorCovarianceBUMP<MODEL>> {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment<MODEL>   Increment_;
  typedef OoBump<MODEL>            OoBump_;
  typedef oops::State<MODEL>       State_;
  typedef ParametersBUMP<MODEL>    Parameters_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceBUMP";}

  ErrorCovarianceBUMP(const Geometry_ &, const oops::Variables &,
                      const eckit::Configuration &, const State_ &, const State_ &);
  virtual ~ErrorCovarianceBUMP();

 private:
  ErrorCovarianceBUMP(const ErrorCovarianceBUMP&);
  ErrorCovarianceBUMP& operator=(const ErrorCovarianceBUMP&);

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<OoBump_> ooBump_;
};

// =============================================================================

template<typename MODEL>
ErrorCovarianceBUMP<MODEL>::ErrorCovarianceBUMP(const Geometry_ & resol,
                                                const oops::Variables & vars,
                                                const eckit::Configuration & conf,
                                                const State_ & xb, const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf), ooBump_()
{
  oops::Log::trace() << "ErrorCovarianceBUMP::ErrorCovarianceBUMP starting" << std::endl;

// Setup timeslots
  std::vector<util::DateTime> timeslots;
  timeslots.push_back(xb.validTime());

// Setup parameters
  Parameters_ param(resol, vars, timeslots, conf);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump_(param.getOoBump()));

  oops::Log::trace() << "ErrorCovarianceBUMP::ErrorCovarianceBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP() {
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovarianceBUMP");
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::~ErrorCovarianceBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doRandomize(Increment_ & dx) const {
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  ooBump_->randomize(dx);
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doMultiply(const Increment_ & dxi,
                                            Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  ooBump_->multiplyNicas(dxi, dxo);
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                                   Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  ooBump_->inverseMultiplyNicas(dxi, dxo);
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceBUMP<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovarianceBUMP<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovarianceBUMP<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCEBUMP_H_
