/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QUENCH_FIELDS_H_
#define QUENCH_FIELDS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace quench {
  class Geometry;

// -----------------------------------------------------------------------------
/// Class to represent a Fields for the  model
class Fields : public util::Printable,
                 public util::Serializable,
                 private util::ObjectCounter<Fields> {
 public:
  static const std::string classname() {return "quench::Fields";}

// Constructors and basic operators
  Fields(const Geometry &, const oops::Variables &, const util::DateTime &);

  void zero();
  void dirac(const eckit::Configuration &);

// ATLAS FieldSet
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}

  const util::DateTime & time() const {return time_;}

/// Serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::shared_ptr<const Geometry> geom_;
  const oops::Variables vars_;
  util::DateTime time_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
};
// -----------------------------------------------------------------------------

}  // namespace quench
#endif  // QUENCH_FIELDS_H_