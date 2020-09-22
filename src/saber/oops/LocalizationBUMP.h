/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_LOCALIZATIONBUMP_H_
#define SABER_OOPS_LOCALIZATIONBUMP_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/LocalizationBase.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "saber/oops/OoBump.h"
#include "saber/oops/ParametersBUMP.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// BUMP localization matrix.

template<typename MODEL> class LocalizationBUMP : public oops::LocalizationBase<MODEL> {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::Increment4D<MODEL>                        Increment4D_;
  typedef OoBump<MODEL>                                   OoBump_;
  typedef ParametersBUMP<MODEL>                           Parameters_;
 public:
  LocalizationBUMP(const Geometry_ &,
                   const eckit::Configuration &);
  ~LocalizationBUMP();

  void multiply(Increment_ &) const override;
  void multiply(Increment4D_ &) const override;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<OoBump_> ooBump_;
  std::vector<util::DateTime> timeslots_;
};

// =============================================================================

template<typename MODEL>
LocalizationBUMP<MODEL>::LocalizationBUMP(const Geometry_ & resol,
                                          const eckit::Configuration & conf)
  : ooBump_()
{
// Setup variables
  const oops::Variables vars(conf, "localization variables");

//  Setup timeslots
  const std::vector<std::string> timeslots_str(conf.getStringVector("timeslots"));
  std::vector<util::DateTime> timeslots;
  for (unsigned jsub = 0; jsub < timeslots_str.size(); ++jsub) {
    const util::DateTime date(timeslots_str[jsub]);
    timeslots.push_back(date);
  }
  oops::Log::info() << "Number of ensemble time-slots:" << timeslots.size() << std::endl;

// Setup parameters
  Parameters_ param(resol, vars, timeslots, conf);

// Transfer OoBump pointer
  ooBump_.reset(new OoBump_(param.getOoBump()));

  oops::Log::trace() << "LocalizationBUMP:LocalizationBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
LocalizationBUMP<MODEL>::~LocalizationBUMP() {
  oops::Log::trace() << "LocalizationBUMP:~LocalizationBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment_ & dx) const {
  oops::Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;
  ooBump_->multiplyNicas(dx);
  oops::Log::trace() << "LocalizationBUMP:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::multiply(Increment4D_ & dx) const {
  oops::Log::trace() << "LocalizationBUMP:multiply starting" << std::endl;
  ooBump_->multiplyNicas(dx);
  oops::Log::trace() << "LocalizationBUMP:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void LocalizationBUMP<MODEL>::print(std::ostream & os) const {
  os << "LocalizationBUMP:print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_LOCALIZATIONBUMP_H_
