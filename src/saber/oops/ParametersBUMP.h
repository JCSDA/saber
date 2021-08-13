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

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/runtime/Exception.h"

#include "eckit/config/Configuration.h"

#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "saber/oops/OoBump.h"

using atlas::Field;
using atlas::FieldSet;
using atlas::array::make_view;

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
                 const oops::Variables &,
                 const eckit::Configuration &,
                 const EnsemblePtr_ ens1 = NULL,
                 const EnsemblePtr_ ens2 = NULL);
  ~ParametersBUMP();

  OoBump_ & getOoBump() {return *ooBump_;}
  void write() const;
  void apply() const;

 private:
  const Geometry_ resol_;
  const oops::Variables inputVars_;
  const oops::Variables activeVars_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<OoBump_> ooBump_;
};

// =============================================================================

template<typename MODEL>
ParametersBUMP<MODEL>::ParametersBUMP(const Geometry_ & resol,
                                      const oops::Variables & inputVars,
                                      const oops::Variables & activeVars,
                                      const eckit::Configuration & conf,
                                      const EnsemblePtr_ ens1,
                                      const EnsemblePtr_ ens2)
  : resol_(resol), inputVars_(inputVars), activeVars_(activeVars), conf_(conf),
    ooBump_() {
  oops::Log::trace() << "ParametersBUMP<MODEL>::ParametersBUMP construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersBUMP");

  // Setup BUMP configuration
  eckit::LocalConfiguration BUMPConf(conf_, "bump");

  // Get ensemble 1 size if ensemble 1 is available
  int ens1_ne = 0;
  if (ens1) ens1_ne = ens1->size();
  if (BUMPConf.has("ensemble")) {
    const eckit::LocalConfiguration ensembleConfig(BUMPConf, "ensemble");
    std::vector<eckit::LocalConfiguration> memberConfig;
    ensembleConfig.get("members", memberConfig);
    ens1_ne = memberConfig.size();
  }
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

  // Read universe size
  oops::Log::info() << "Read universe radius" << std::endl;
  std::unique_ptr<atlas::FieldSet> universe_rad(new atlas::FieldSet());
  if (conf_.has("universe radius")) {
    // Get configuration
    const eckit::LocalConfiguration universeRadiusConf(conf_, "universe radius");

    // Setup increment
    util::DateTime time(1977, 5, 25, 0, 0, 0);
    Increment_ dx(resol_, activeVars_, time);
    dx.read(universeRadiusConf);

    // Get ATLAS fieldset
    dx.toAtlas(universe_rad.get());
  }

  // Create BUMP
  oops::Log::info() << "Create BUMP" << std::endl;
  ooBump_.reset(new OoBump_(resol, activeVars_, BUMPConf, universe_rad.get()));

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
      // Get date
      const util::DateTime date(inputConf.getString("date"));

      // Setup increment
      Increment_ dx(resol_, activeVars_, date);
      dx.read(inputConf);

      // Set parameter to BUMP
      std::string param = inputConf.getString("parameter");
      ooBump_->setParameter(param, dx);
      oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                        << std::setprecision(3) << dx.norm() << std::endl;
    }
  }

  // Load ensemble members sequentially
  if (BUMPConf.has("ensemble")) {
    // Get ensemble and members configurations
    const eckit::LocalConfiguration ensembleConfig(BUMPConf, "ensemble");
    std::vector<eckit::LocalConfiguration> memberConfig;
    ensembleConfig.get("members", memberConfig);

    // Check what needs to be updated
    int update_vbal_cov = BUMPConf.getInt("update_vbal_cov", 0);
    int update_var = BUMPConf.getInt("update_var", 0);
    int update_mom = BUMPConf.getInt("update_mom", 0);

    // Loop over all ensemble members
    for (int ie = 0; ie < ens1_ne; ++ie) {
      // Get date
      const util::DateTime date(memberConfig[ie].getString("date"));

      // Define increment
      Increment_ incr(resol, activeVars_, date);

      // Read member
      oops::Log::info() <<
      "-------------------------------------------------------------------" << std::endl;
      oops::Log::info() << "--- Load member " << ie+1 << " / " << ens1_ne << std::endl;
      incr.read(memberConfig[ie]);

      if (update_vbal_cov == 1) {
        // Update vertical covariance
        ooBump_->updateVbalCov(incr, ie);
      }
      if (update_var == 1) {
        // Update variance
        ooBump_->updateVar(incr, ie);
      }
      if (update_mom == 1) {
        // Update moments
        ooBump_->updateMom(incr, ie);
      }
    }
  }

  // Estimate parameters
  ooBump_->runDrivers();

  // Partial deallocation
  ooBump_->partialDealloc();

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
  if (outputConfs.size() > 0) {
    for (const auto & outputConf : outputConfs) {
      // Get date
      const util::DateTime date(outputConf.getString("date"));

      // Setup increment
      Increment_ dx(resol_, activeVars_, date);

      // Set increment to zero
      dx.zero();

      // Get parameter from BUMP
      std::string param = outputConf.getString("parameter");
      ooBump_->getParameter(param, dx);

      // Write parameter
      dx.write(outputConf);
      oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                        << std::setprecision(3) << dx.norm() << std::endl;
    }
  } else {
    oops::Log::test() << "No output configuration" << std::endl;
  }

  oops::Log::trace() << "ParametersBUMP::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersBUMP<MODEL>::apply() const {
  oops::Log::trace() << "ParametersBUMP::apply starting" << std::endl;
  util::Timer timer(classname(), "apply");

  // Write parameters
  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::info() << "--- Apply operators" << std::endl;

  std::vector<eckit::LocalConfiguration> appConfs;
  conf_.get("operators application", appConfs);
  if (appConfs.size() > 0) {
    for (const auto & appConf : appConfs) {
      // Get date
      const util::DateTime date(appConf.getString("date"));

      // Setup increments
      Increment_ dxi(resol_, inputVars_, date);
      Increment_ dxo(resol_, inputVars_, date);

      // Read input file
      eckit::LocalConfiguration inputConf(appConf, "input");
      oops::Log::info() << "       - Input file: " << inputConf << std::endl;
      dxi.read(inputConf);

      // Apply BUMP operator
      std::vector<std::string> bumpOperators;
      appConf.get("bump operators", bumpOperators);
      for (const auto & bumpOperator : bumpOperators) {
        oops::Log::info() << "         Apply " << bumpOperator << std::endl;
        if (bumpOperator == "multiplyVbal") ooBump_->multiplyVbal(dxi, dxo);
        if (bumpOperator == "multiplyVbalInv") ooBump_->multiplyVbalInv(dxi, dxo);
        if (bumpOperator == "multiplyVbalAd") ooBump_->multiplyVbalAd(dxi, dxo);
        if (bumpOperator == "multiplyVbalInvAd") ooBump_->multiplyVbalInvAd(dxi, dxo);
        if (bumpOperator == "multiplyStdDev") ooBump_->multiplyStdDev(dxi, dxo);
        if (bumpOperator == "multiplyStdDevInv") ooBump_->multiplyStdDevInv(dxi, dxo);
        if (bumpOperator == "multiplyNicas") ooBump_->multiplyNicas(dxi, dxo);
        if (bumpOperator == "inverseMultiplyNicas") ooBump_->inverseMultiplyNicas(dxi, dxo);
// TODO(Benjamin): psichi_to_uv ?
        dxi = dxo;
      }

      // Write file
      eckit::LocalConfiguration outputConf(appConf, "output");
      oops::Log::info() << "         Output file: " << outputConf << std::endl;
      dxo.write(outputConf);
    }
  }

  oops::Log::info() <<
  "-------------------------------------------------------------------" << std::endl;
  oops::Log::trace() << "ParametersBUMP::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_PARAMETERSBUMP_H_
