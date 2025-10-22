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

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {
  class Fields;

// -----------------------------------------------------------------------------
/// FieldsIOBase class

class FieldsIOBase : private eckit::NonCopyable {
 public:
  static const std::string classname()
    {return "quench::FieldsIOBase";}

  // Constructor
  explicit FieldsIOBase(const std::string & ioFormat)
    : ioFormat_(ioFormat) {}

  // Read
  virtual void read(const oops::Variables &,
                    const eckit::Configuration &,
                    Fields &) const
    {throw eckit::Exception("read not implemented for this format", Here());}

  // Write
  virtual void write(const eckit::Configuration &,
                     const Fields &) const
    {throw eckit::Exception("read not implemented for this format", Here());}

 protected:
  const std::string ioFormat_;
};

// -----------------------------------------------------------------------------
/// FieldsIOFactory class

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
/// FieldsIOMaker class

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
