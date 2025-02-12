/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Increment.h"

#include <vector>

#include "atlas/field.h"

#include "oops/util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const oops::Variables & vars,
                     const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)) {
  oops::Log::trace() << classname() << "::Increment starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << classname() << "::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment::Increment(const Geometry & resol,
                     const Increment & other)
  : fields_(new Fields(*other.fields_, resol)) {
  oops::Log::trace() << classname() << "::Increment starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << classname() << "::Increment done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment::Increment(const Increment & other,
                     const bool copy)
  : fields_(new Fields(*other.fields_, copy)) {
  oops::Log::trace() << classname() << "::Increment" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::diff(const State & x1,
                     const State & x2) {
  oops::Log::trace() << classname() << "::diff starting" << std::endl;

  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());

  oops::Log::trace() << classname() << "::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator=(const Increment & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  fields_.reset(new Fields(*rhs.fields_));

  oops::Log::trace() << classname() << "::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator+=(const Increment & dx) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator-=(const Increment & dx) {
  oops::Log::trace() << classname() << "::operator-= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;

  oops::Log::trace() << classname() << "::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

Increment & Increment::operator*=(const double & zz) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  *fields_ *= zz;

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void Increment::zero(const util::DateTime & vt) {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  fields_->zero();
  fields_->time() = vt;

  oops::Log::trace() << classname() << "::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::axpy(const double & zz,
                     const Increment & dx,
                     const bool check) {
  oops::Log::trace() << classname() << "::axpy starting" << std::endl;

  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);

  oops::Log::trace() << classname() << "::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

void Increment::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << std::endl << "- Valid time: " << this->validTime();
  os << *fields_;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
