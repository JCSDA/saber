/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_PARAMETERSBUMP_H_
#define SABER_OOPS_PARAMETERSBUMP_H_

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "saber/oops/OoBump.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// BUMP parameters

template<typename MODEL>
class ParametersBUMP {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef OoBump<MODEL>                                   OoBump_;
  typedef oops::State<MODEL>                              State_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::ParametersBUMP";}
  ParametersBUMP(const Geometry_ &,
                 const oops::Variables &,
                 const util::DateTime &,
                 const eckit::Configuration &,
                 const EnsemblePtr_ ens1 = NULL,
                 const EnsemblePtr_ ens2 = NULL);
  ~ParametersBUMP();

  OoBump_ & getOoBump() {return *ooBump_;}
  void write() const;

 private:
  const Geometry_ resol_;
  const oops::Variables vars_;
  util::DateTime time_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<OoBump_> ooBump_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const Geometry_ & resol,
                                      const oops::Variables & vars,
                                      const util::DateTime & time,
                                      const eckit::Configuration & conf,
                                      const EnsemblePtr_ ens1,
                                      const EnsemblePtr_ ens2)
  : resol_(resol), vars_(vars), time_(time), conf_(conf), ooBump_()
{
  oops::Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

  // Setup BUMP configuration
  eckit::LocalConfiguration BUMPConf(conf_, "bump");

  // Get ensemble 1 size if ensemble 1 is available
  int ens1_ne = 0;
  if (ens1) ens1_ne = ens1->size();
  BUMPConf.set("ens1_ne", ens1_ne);
  if (!BUMPConf.has("ens1_nsub")) BUMPConf.set("ens1_nsub", 1);

  // Get ensemble 2 size if ensemble 2 is available
  int ens2_ne = 0;
  if (ens2) ens2_ne = ens2->size();
  BUMPConf.set("ens2_ne", ens2_ne);
  if (!BUMPConf.has("ens2_nsub")) BUMPConf.set("ens2_nsub", 1);

  // Get missing value
  const double msvalr = util::missingValue(msvalr);
  BUMPConf.set("msvalr", msvalr);

  // Create BUMP
  oops::Log::info() << "Create BUMP" << std::endl;
  ooBump_.reset(new OoBump_(resol, vars, time_, BUMPConf));

  // Add members of ensemble 1
  if (ens1) {
    oops::Log::info() << "--- Add members of ensemble 1" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      oops::Log::info() << "      Member " << ie+1 << " / " << ens1_ne << std::endl;
      ooBump_->addMember((*ens1)[ie].atlas(), ie, 1);
    }
  }

  // Add members of ensemble 2
  if (ens2) {
    oops::Log::info() << "--- Add members of ensemble 2" << std::endl;
    for (int ie = 0; ie < ens2_ne; ++ie) {
      oops::Log::info() << "      Member " << ie+1 << " / " << ens2_ne << std::endl;
      ooBump_->addMember((*ens2)[ie].atlas(), ie, 2);
    }
  }

  // Read data from files
  oops::Log::info() << "    Read data from files" << std::endl;
  if (conf_.has("input")) {
    // Set BUMP input parameters
    std::vector<eckit::LocalConfiguration> inputConfs;
    conf_.get("input", inputConfs);

    for (const auto & inputConf : inputConfs) {
      // Read parameter for the specified time
      const util::DateTime date(inputConf.getString("date"));

      // Setup increment
      Increment_ dx(resol_, vars_, time_);
      dx.read(inputConf);

      // Set parameter to BUMP
      std::string param = inputConf.getString("parameter");
      ooBump_->setParameter(param, dx);
    }
  }

  // Estimate parameters
  ooBump_->runDrivers();

  oops::Log::trace() << "ParametersBUMP:ParametersBUMP constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ParametersBUMP<MODEL>::~ParametersBUMP() {
  oops::Log::trace() << "ParametersBUMP<MODEL>::~ParametersBUMP destruction starting" << std::endl;
  util::Timer timer(classname(), "~ParametersBUMP");
  oops::Log::trace() << "ParametersBUMP:~ParametersBUMP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersBUMP<MODEL>::write() const {
  oops::Log::trace() << "ParametersBUMP::write starting" << std::endl;
  util::Timer timer(classname(), "write");

// Write parameters
  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::info() << "--- Write parameters" << std::endl;

  std::vector<eckit::LocalConfiguration> outputConfs;
  conf_.get("output", outputConfs);
  for (const auto & outputConf : outputConfs) {
  // Setup dummy increment
    Increment_ dx(resol_, vars_, time_);
    dx.zero();

  // Get parameter from BUMP
    std::string param = outputConf.getString("parameter");
    ooBump_->getParameter(param, dx);

  // Write parameter for the specified time
    const util::DateTime date(outputConf.getString("date"));
    dx.write(outputConf);
    oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                      << std::setprecision(3) << dx.norm() << std::endl;
  }
  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::trace() << "ParametersBUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_PARAMETERSBUMP_H_
