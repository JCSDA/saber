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
void Increment::zero() {
  fields_->zero();
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
  os << std::endl << "Valid time:" << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------

}  // namespace quench
