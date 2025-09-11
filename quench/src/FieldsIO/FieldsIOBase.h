/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/base/Variables.h"

namespace quench {
  class Geometry;

// -----------------------------------------------------------------------------

class FieldsIOBase : private eckit::NonCopyable {
 public:
  static const std::string classname()
    {return "quench::FieldsIOBase";}

  // Constructor/destructor
  explicit FieldsIOBase(const std::string & ioFormat)
    : ioFormat_(ioFormat) {}

  // Read
  virtual void read(const Geometry &,
                    const oops::Variables &,
                    const eckit::Configuration &,
                    atlas::FieldSet &) const
    {throw eckit::Exception("read not implemented for this format", Here());}

  // Write
  virtual void write(const Geometry &,
                     const eckit::Configuration &,
                     const atlas::FieldSet &) const
    {throw eckit::Exception("read not implemented for this format", Here());}

 protected:
  const std::string ioFormat_;
};

// -----------------------------------------------------------------------------

class FieldsIOFactory;

// -----------------------------------------------------------------------------

class FieldsIOFactory {
 public:
  static const std::string classname()
    {return "quench::FieldsIOFactory";}

  static std::unique_ptr<FieldsIOBase> create(const std::string &);

  virtual ~FieldsIOFactory() = default;

 protected:
  explicit FieldsIOFactory(const std::string &);

 private:
  virtual std::unique_ptr<FieldsIOBase> make(const std::string &) = 0;

  static std::map < std::string, FieldsIOFactory * > & getMakers() {
    static std::map < std::string, FieldsIOFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class FieldsIOMaker : public FieldsIOFactory {
  std::unique_ptr<FieldsIOBase> make(const std::string & ioFormat) override {
    return std::make_unique<T>(ioFormat);
  }

 public:
  explicit FieldsIOMaker(const std::string & ioFormat) : FieldsIOFactory(ioFormat) {}
};

// -----------------------------------------------------------------------------

}  // namespace quench
