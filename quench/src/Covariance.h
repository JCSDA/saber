/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/Geometry.h"

namespace quench {
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for quench model.

class Covariance : public util::Printable,
                   private boost::noncopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname() {return "quench::Covariance";}

  Covariance(const Geometry &, const oops::Variables &,
             const eckit::Configuration &, const State &, const State &)
    {throw eckit::NotImplemented(Here());}

  void multiply(const Increment &, Increment &) const
    {throw eckit::NotImplemented(Here());}
  void inverseMultiply(const Increment &, Increment &) const
    {throw eckit::NotImplemented(Here());}
  void randomize(Increment &) const
    {throw eckit::NotImplemented(Here());}

 private:
  void print(std::ostream & os) const {os << "Covariance";}
};
// -----------------------------------------------------------------------------

}  // namespace quench
