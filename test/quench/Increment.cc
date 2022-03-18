/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quench/Increment.h"

#include <vector>

#include "atlas/field.h"

#include "oops/util/Logger.h"

#include "quench/Fields.h"

namespace quench {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & resol, const oops::Variables & vars,
                     const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt))
{
  fields_->zero();
  oops::Log::trace() << "Increment constructed." << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  : fields_(new Fields(*other.fields_, copy))
{
  oops::Log::trace() << "Increment copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());
}
// -----------------------------------------------------------------------------
Increment & Increment::operator=(const Increment & rhs) {
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator*=(const double & zz) {
  *fields_ *= zz;
  return *this;
}

// -----------------------------------------------------------------------------
void Increment::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void Increment::axpy(const double & zz, const Increment & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);
}
// -----------------------------------------------------------------------------
void Increment::accumul(const double & zz, const State & xx) {
  fields_->axpy(zz, xx.fields());
}
// -----------------------------------------------------------------------------
void Increment::schur_product_with(const Increment & dx) {
  fields_->schur_product_with(*dx.fields_);
}
// -----------------------------------------------------------------------------
double Increment::dot_product_with(const Increment & other) const {
  return fields_->dot_product_with(*other.fields_);
}
// -----------------------------------------------------------------------------
void Increment::random() {
  fields_->random();
}
// -----------------------------------------------------------------------------
void Increment::dirac(const eckit::Configuration & config) {
  fields_->dirac(config);
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void Increment::setAtlas(atlas::FieldSet * afieldset) const {
  fields_->setAtlas(afieldset);
}
// -----------------------------------------------------------------------------
void Increment::toAtlas(atlas::FieldSet * afieldset) const {
  fields_->toAtlas(afieldset);
}
// -----------------------------------------------------------------------------
void Increment::fromAtlas(atlas::FieldSet * afieldset) {
  fields_->fromAtlas(afieldset);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void Increment::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void Increment::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t Increment::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void Increment::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void Increment::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void Increment::print(std::ostream & os) const {
  os << std::endl << "Valid time:" << this->validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------

}  // namespace quench
