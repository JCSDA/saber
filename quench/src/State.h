/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/Fields.h"
#include "src/Increment.h"

namespace quench {
  class Geometry;
  class Increment;

/// State class
// -----------------------------------------------------------------------------
class State : public util::Printable,
                private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "quench::State";}

/// Constructor, destructor
  State(const Geometry &, const oops::Variables &, const util::DateTime &);
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const oops::Variables &, const State &);
  State(const State &);

/// Assignment
  State & operator=(const State &);

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {return fields_->geometry();}
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const;;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

/// ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

/// Other
  void zero();
  void accumul(const double &, const State &);

 private:
  void print(std::ostream &) const;
  std::unique_ptr<Fields> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
