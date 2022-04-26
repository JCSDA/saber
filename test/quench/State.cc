/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "quench/State.h"

#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "quench/Fields.h"

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
  oops::Variables vars = oops::Variables(file, "state variables");
  fields_.reset(new Fields(resol, vars, util::DateTime()));
  fields_->zero();
  const util::DateTime vt(file.getString("date"));
  fields_->time() = vt;
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry &, const State &) : fields_() {
  throw eckit::NotImplemented("State constructor", Here());
}
// -----------------------------------------------------------------------------
State::State(const State & other)
  : fields_(new Fields(*other.fields_))
{
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
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
