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

#include "oops/assimilation/Increment4D.h"
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
  typedef oops::Increment4D<MODEL>                        Increment4D_;
  typedef OoBump<MODEL>                                   OoBump_;
  typedef oops::State<MODEL>                              State_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  static const std::string classname() {return "oops::ParametersBUMP";}
  ParametersBUMP(const Geometry_ &,
                 const oops::Variables &,
                 const std::vector<util::DateTime> &,
                 const eckit::Configuration &,
                 const EnsemblePtr_ ens = NULL,
                 const EnsemblePtr_ pseudo_ens = NULL);
  ~ParametersBUMP();

  OoBump_ & getOoBump() {return *ooBump_;}
  void write() const;

 private:
  const Geometry_ resol_;
  const oops::Variables vars_;
  std::vector<util::DateTime> timeslots_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<OoBump_> ooBump_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const Geometry_ & resol,
                                      const oops::Variables & vars,
                                      const std::vector<util::DateTime> & timeslots,
                                      const eckit::Configuration & conf,
                                      const EnsemblePtr_ ens,
                                      const EnsemblePtr_ pseudo_ens)
  : resol_(resol), vars_(vars), timeslots_(timeslots), conf_(conf), ooBump_()
{
  oops::Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

  // Setup BUMP configuration
  eckit::LocalConfiguration BUMPConf(conf_, "bump");

  // Setup members release
  int release_members = 0;
  if (BUMPConf.has("release_members")) release_members = BUMPConf.getInt("release_members");

  // Get ensemble size if ensemble is available
  int ens1_ne = 0;
  if (ens) ens1_ne = ens->size();
  BUMPConf.set("ens1_ne", ens1_ne);
  BUMPConf.set("ens1_nsub", 1);

  // Get pseudo-ensemble size if pseudo-ensemble is available
  int ens2_ne = 0;
  if (pseudo_ens) ens2_ne = pseudo_ens->size();
  BUMPConf.set("ens2_ne", ens2_ne);
  BUMPConf.set("ens2_nsub", 1);

  // Get missing value
  const double msvalr = util::missingValue(msvalr);
  BUMPConf.set("msvalr", msvalr);

  // Create BUMP
  oops::Log::info() << "Create BUMP" << std::endl;
  ooBump_.reset(new OoBump_(resol, vars, timeslots_, BUMPConf));

// Transfer/copy ensemble members to BUMP
  if (release_members == 1) {
    oops::Log::info() << "Transfer ensemble members to BUMP" << std::endl;
  } else {
    oops::Log::info() << "Copy ensemble members to BUMP" << std::endl;
  }
  for (int ie = 0; ie < ens1_ne; ++ie) {
    oops::Log::info() << "   Member " << ie+1 << " / " << ens1_ne << std::endl;;

  // Setup increment
    Increment4D_ dx(resol_, vars_, timeslots_);

  // Copy member
    if (release_members == 1) {
      dx = (*ens)[0];
    } else {
      dx = (*ens)[ie];
    }

  // Copy data to BUMP
    ooBump_->addMember(dx, ie);

    if (release_members == 1) {
    // Release ensemble member
      ens->releaseMember();
    }
  }

// Transfer/copy pseudo-ensemble members to BUMP
  if (release_members == 1) {
    oops::Log::info() << "Transfer pseudo-ensemble members to BUMP" << std::endl;
  } else {
    oops::Log::info() << "Copy pseudo-ensemble members to BUMP" << std::endl;
  }
  for (int ie = 0; ie < ens2_ne; ++ie) {
    oops::Log::info() << "   Member " << ie+1 << " / " << ens2_ne << std::endl;

  // Setup increment
    Increment4D_ dx(resol_, vars_, timeslots_);

  // Copy member
    if (release_members == 1) {
      dx = (*pseudo_ens)[0];
    } else {
      dx = (*pseudo_ens)[ie];
    }

  // Copy data to BUMP
    ooBump_->addPseudoMember(dx, ie);

    if (release_members == 1) {
    // Release ensemble member
      pseudo_ens->releaseMember();
    }
  }

// Read data from files
  oops::Log::info() << "Read data from files" << std::endl;
  if (conf_.has("input")) {
  // Set BUMP input parameters
    std::vector<eckit::LocalConfiguration> inputConfs;
    conf_.get("input", inputConfs);

    for (const auto & inputConf : inputConfs) {
    // Read parameter for the specified timeslots
      const util::DateTime date(inputConf.getString("date"));
      bool found = false;

    // Setup increment
      Increment4D_ dx(resol_, vars_, timeslots_);
      dx.zero();

      for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
        if (date == timeslots_[jsub]) {
          found = true;
          dx[dx.first()+jsub].read(inputConf);
        }
      }
      ASSERT(found);

    // Set parameter to BUMP
      std::string param = inputConf.getString("parameter");
      ooBump_->setParameter(param, dx);
    }
  }

// Estimate parameters
  ooBump_->runDrivers();

  if (release_members == 1) {
  // Transfer ensemble members from BUMP
    oops::Log::info() << "Transfer ensemble members from BUMP" << std::endl;
    for (int ie = 0; ie < ens1_ne; ++ie) {
      oops::Log::info() << "   Member " << ie+1 << " / " << ens1_ne << std::endl;


    // Setup dummy increment
      Increment4D_ dx(resol_, vars_, timeslots_);

    // Copy data from BUMP
      ooBump_->removeMember(dx, ie);

    // Reset ensemble member
      ens->resetMember(dx);
    }

  // Transfer pseudo-ensemble members from BUMP
    oops::Log::info() << "Transfer pseudo-ensemble members from BUMP" << std::endl;
    for (int ie = 0; ie < ens2_ne; ++ie) {
      oops::Log::info() << "   Member " << ie+1 << " / " << ens2_ne << std::endl;

    // Setup dummy increment
      Increment4D_ dx(resol_, vars_, timeslots_);

    // Copy data from BUMP
      ooBump_->removePseudoMember(dx, ie);

    // Reset pseudo-ensemble member
      pseudo_ens->resetMember(dx);
    }
  }

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
  std::vector<eckit::LocalConfiguration> outputConfs;
  conf_.get("output", outputConfs);
  for (const auto & outputConf : outputConfs) {
  // Setup dummy increment
    Increment4D_ dx(resol_, vars_, timeslots_);
    dx.zero();

  // Get parameter from BUMP
    std::string param = outputConf.getString("parameter");
    ooBump_->getParameter(param, dx);

  // Write parameter for the specified timeslots
    const util::DateTime date(outputConf.getString("date"));
    bool found = false;
    for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
      int isub = jsub+dx.first();
      if (date == timeslots_[jsub]) {
        found = true;
        dx[isub].write(outputConf);
        oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                    << std::setprecision(3) << dx[isub].norm() << std::endl;
      }
    }
    ASSERT(found);
  }
  oops::Log::trace() << "ParametersBUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_PARAMETERSBUMP_H_
