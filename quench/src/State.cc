/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/State.h"

#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "src/Fields.h"
#include "src/Increment.h"

namespace quench {

// -----------------------------------------------------------------------------

State::State(const Geometry & resol,
             const oops::Variables & vars,
             const util::DateTime & vt)
  : fields_(new Fields(resol, vars, vt)) {
  oops::Log::trace() << classname() << "::State starting" << std::endl;

  fields_->zero();

  oops::Log::trace() << classname() << "::State done" << std::endl;
}

// -----------------------------------------------------------------------------

State::State(const Geometry & resol,
             const eckit::Configuration & file)
  : fields_() {
  oops::Log::trace() << classname() << "::State starting" << std::endl;

  const std::vector<std::string> varNames = file.has("state variables") ?
    file.getStringVector("state variables") : file.getStringVector("variables");
  const oops::Variables vars(varNames);
  fields_.reset(new Fields(resol, vars, util::DateTime()));
  if (file.has("filepath")) {
    oops::Log::info() << "Info     : Create state from file" << std::endl;
    fields_->read(file);
  } else if (file.has("constant value")) {
    oops::Log::info() << "Info     : Create state with a constant value" << std::endl;
    fields_->constantValue(file.getDouble("constant value"));
  } else if (file.has("constant profile")) {
    oops::Log::info() << "Info     : Create state with a constant profile" << std::endl;
    fields_->constantValue(file.getDoubleVector("constant profile"));
  } else if (file.has("constant group-specific value")) {
    oops::Log::info() << "Info     : Create state with a constant group-specific value"
                      << std::endl;
    fields_->constantValue(file);
  } else {
    oops::Log::info() << "Info     : Create empty state" << std::endl;
    fields_->zero();
  }
  const util::DateTime vt(file.getString("date"));
  fields_->time() = vt;

  oops::Log::trace() << classname() << "::State done" << std::endl;
}

// -----------------------------------------------------------------------------

State & State::operator=(const State & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  fields_.reset(new Fields(*rhs.fields_));

  oops::Log::trace() << classname() << "::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

State & State::operator+=(const Increment & dx) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  *fields_+=dx.fields();

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void State::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << std::endl << "- Valid time: " << this->validTime();
  os << *fields_;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quench
