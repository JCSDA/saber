/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_LOCALIZATIONID_H_
#define SABER_OOPS_LOCALIZATIONID_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/LocalizationBase.h"
#include "oops/util/Logger.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// Identity localization matrix.

template<typename MODEL> class LocalizationID : public oops::LocalizationBase<MODEL> {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;

 public:
  LocalizationID(const Geometry_ &,
                 const util::DateTime & time,
                 const eckit::Configuration &);
  ~LocalizationID();

  void doRandomize(Increment_ &) const override;
  void doMultiply(Increment_ &) const override;

 private:
  void print(std::ostream &) const override;
};

// =============================================================================

template<typename MODEL>
LocalizationID<MODEL>::LocalizationID(const Geometry_ & resol,
                                      const util::DateTime & time,
                                      const eckit::Configuration & conf) {
  oops::Log::trace() << "LocalizationID:LocalizationID constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationID<MODEL>::~LocalizationID() {
  oops::Log::trace() << "LocalizationID:~LocalizationID destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::doRandomize(Increment_ & dx) const {
  dx.random();
  oops::Log::trace() << "LocalizationID:doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::doMultiply(Increment_ & dx) const {
  oops::Log::trace() << "LocalizationID:doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationID<MODEL>::print(std::ostream & os) const {
  os << "LocalizationID:print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_LOCALIZATIONID_H_
