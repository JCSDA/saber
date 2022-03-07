/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_INCREMENT_H_
#define QUENCH_INCREMENT_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "quench/Fields.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Geometry;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class Increment : public util::Printable,
                    public util::Serializable,
                    private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "quench::Increment";}

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &, const util::DateTime &);

/// Basic operators
  void zero();
  void dirac(const eckit::Configuration &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}

/// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const {return fields_->geometry();}

/// Other
  const oops::Variables & variables() const {return fields_->variables();}

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

/// Data
 private:
  void print(std::ostream &) const override;
  std::unique_ptr<Fields> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace quench

#endif  // QUENCH_INCREMENT_H_
