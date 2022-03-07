/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_STATE_H_
#define QUENCH_STATE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quench/Fields.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Geometry;

/// State class
// -----------------------------------------------------------------------------
class State : public util::Printable,
                private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "quench::State";}

/// Constructor, destructor
  State(const Geometry &, const oops::Variables &, const util::DateTime &);

/// I/O and diagnostics
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {return fields_->geometry();}
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const;;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// Other
  void zero();

 private:
  void print(std::ostream &) const;
  std::unique_ptr<Fields> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace quench

#endif  // QUENCH_STATE_H_
