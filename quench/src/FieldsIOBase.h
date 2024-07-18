/*
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"
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
  explicit FieldsIOBase(const std::string &)
    {}

  // Read
  virtual void read(const Geometry &,
                    const oops::Variables &,
                    const eckit::Configuration &,
                    atlas::FieldSet &) const = 0;

  // Write
  virtual void write(const Geometry &,
                     const eckit::Configuration &,
                     const atlas::FieldSet &) const = 0;
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
