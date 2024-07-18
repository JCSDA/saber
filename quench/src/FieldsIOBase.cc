/*
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 */

#include "src/FieldsIOBase.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace quench {

// -----------------------------------------------------------------------------

FieldsIOFactory::FieldsIOFactory(const std::string & ioFormat) {
  if (getMakers().find(ioFormat) != getMakers().end()) {
    oops::Log::error() << ioFormat << " already registered in quench::FieldsIOFactory."
      << std::endl;
    throw eckit::Exception("Element already registered in quench::FieldsIOFactory.",
      Here());
  }
  getMakers()[ioFormat] = this;
}

// -----------------------------------------------------------------------------

std::unique_ptr<FieldsIOBase> FieldsIOFactory::create(const std::string & ioFormat) {
  oops::Log::trace() << classname() << "::create starting" << std::endl;

  typename std::map<std::string, FieldsIOFactory*>::iterator jsb = getMakers().find(ioFormat);
  if (jsb == getMakers().end()) {
    oops::Log::error() << ioFormat << " does not exist in quench::FieldsIOFactory." << std::endl;
    throw eckit::UserError("Element does not exist in quench::FieldsIOFactory.", Here());
  }
  std::unique_ptr<FieldsIOBase> ptr = jsb->second->make(ioFormat);

  oops::Log::trace() << classname() << "::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace quench
