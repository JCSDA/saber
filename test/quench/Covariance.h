/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_COVARIANCE_H_
#define QUENCH_COVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quench/CovarianceParams.h"
#include "quench/Geometry.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace quench {
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for quench model.

class Covariance : public util::Printable,
                   private boost::noncopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  typedef CovarianceParams Parameters_;
  static const std::string classname() {return "quench::Covariance";}

  Covariance(const Geometry &, const oops::Variables &,
             const Parameters_ &, const State &, const State &) {}

  void multiply(const Increment &, Increment &) const {}
  void inverseMultiply(const Increment &, Increment &) const {}
  void randomize(Increment &) const {}

 private:
  void print(std::ostream &) const {}
};
// -----------------------------------------------------------------------------

}  // namespace quench
#endif  // QUENCH_COVARIANCE_H_
