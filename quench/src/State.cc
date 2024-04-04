/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "src/State.h"

#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "src/Fields.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const oops::Variables & vars,
                 const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt))
{
  fields_->zero();
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const eckit::Configuration & file)
  : fields_()
{
  const oops::Variables vars(file, "state variables");
  fields_.reset(new Fields(resol, vars, util::DateTime()));
  if (file.has("filepath")) {
    oops::Log::info() << "Info     : Create state from file" << std::endl;
    fields_->read(file);
  } else {
    oops::Log::info() << "Info     : Create empty state" << std::endl;
    if (file.has("constant value")) {
      fields_->constantValue(file.getDouble("constant value"));
    } else if (file.has("constant group-specific value")) {
      fields_->constantValue(file);
    } else {
      fields_->zero();
    }
  }
  const util::DateTime vt(file.getString("date"));
  fields_->time() = vt;
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol, const State & other)
  : fields_(new Fields(*other.fields_, resol))
{
  ASSERT(fields_);
  oops::Log::trace() << "State::State created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const State & other)
  : fields_(new Fields(*other.fields_))
{
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  *fields_+=dx.fields();
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void State::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void State::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t State::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void State::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void State::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << this->validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet accessor
// -----------------------------------------------------------------------------
void State::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
// -----------------------------------------------------------------------------
void State::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void State::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void State::accumul(const double & zz, const State & xx) {
  fields_->axpy(zz, *xx.fields_);
}
// -----------------------------------------------------------------------------

}  // namespace quench
