/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
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

#include "src/Increment.h"

namespace quench {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// Covariance class

class Covariance : public util::Printable,
                   private boost::noncopyable,
                   private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname()
    {return "quench::Covariance";}

  // Constructor/destructor
  Covariance(const Geometry &,
             const oops::Variables &,
             const eckit::Configuration &,
             const State &,
             const State &)
    {}
  ~Covariance()
    {}

  // Multiply and inverse multiply (identity)
  void multiply(const Increment & dxi,
                Increment & dxo) const
    {dxo = dxi;}
  void inverseMultiply(const Increment & dxi,
                       Increment & dxo) const
    {dxo = dxi;}

  // Randomization
  void randomize(Increment & dxo) const
    {dxo.random();}

 private:
  // Print
  void print(std::ostream & os) const
    {os << "Covariance";}
};

// -----------------------------------------------------------------------------

}  // namespace quench
