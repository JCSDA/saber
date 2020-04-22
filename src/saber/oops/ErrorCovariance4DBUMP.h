/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCE4DBUMP_H_
#define SABER_OOPS_ERRORCOVARIANCE4DBUMP_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
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

/// Model space error 4D covariance with BUMP

template <typename MODEL>
class ErrorCovariance4DBUMP : public oops::ModelSpaceCovariance4DBase<MODEL>,
                              public util::Printable,
                              private util::ObjectCounter<ErrorCovariance4DBUMP<MODEL> >,
                              private boost::noncopyable {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment4D<MODEL> Increment4D_;
  typedef OoBump<MODEL>            OoBump_;
  typedef oops::State<MODEL>       State_;
  typedef oops::State4D<MODEL>     State4D_;
  typedef ParametersBUMP<MODEL>    Parameters_;

 public:
  static const std::string classname() {return "saber::ErrorCovariance4DBUMP";}

  ErrorCovariance4DBUMP(const Geometry_ &, const oops::Variables &,
                        const eckit::Configuration &, const State4D_ &, const State4D_ &);
  virtual ~ErrorCovariance4DBUMP();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<OoBump_> ooBump_;
  std::vector<util::DateTime> timeslots_;
};

// =============================================================================

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::ErrorCovariance4DBUMP(const Geometry_ & resol,
                                                    const oops::Variables & vars,
                                                    const eckit::Configuration & conf,
                                                    const State4D_ & xb, const State4D_ & fg)
  : oops::ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, conf), ooBump_(),
    timeslots_(xb.validTimes())
{
  oops::Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP starting" << std::endl;

// Setup parameters
  Parameters_ param(resol, vars, timeslots_, conf);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump_(param.getOoBump()));

  oops::Log::trace() << "ErrorCovariance4DBUMP::ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP() {
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP starting"
                     << std::endl;
  util::Timer timer(classname(), "~ErrorCovariance4DBUMP");
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::~ErrorCovariance4DBUMP done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doRandomize(Increment4D_ & dx) const {
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  ooBump_->randomize(dx);
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doMultiply(const Increment4D_ & dxi,
                                              Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  ooBump_->multiplyNicas(dxi, dxo);
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::doInverseMultiply(const Increment4D_ & dxi,
                                                     Increment4D_ & dxo) const {
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  ooBump_->inverseMultiplyNicas(dxi, dxo);
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovariance4DBUMP<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovariance4DBUMP<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovariance4DBUMP<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCE4DBUMP_H_
