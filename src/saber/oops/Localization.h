/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/generic/LocalizationBase.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "saber/oops/ReadInput.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {

// -----------------------------------------------------------------------------

template<typename MODEL>
class Localization : public oops::LocalizationBase<MODEL> {
  typedef oops::Geometry<MODEL>               Geometry_;
  typedef oops::Increment<MODEL>              Increment_;
  typedef oops::State<MODEL>                  State_;

 public:
  Localization(const Geometry_ &,
               const oops::Variables &,
               const eckit::Configuration &);
  ~Localization();

  void randomize(Increment_ &) const override;
  void multiply(Increment_ &) const override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<SaberCentralBlockBase> saberCentralBlock_;
};

// =============================================================================

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geom,
                                  const oops::Variables & incVars,
                                  const eckit::Configuration & conf)
  : saberCentralBlock_()
{
  oops::Log::trace() << "Localization::Localization starting" << std::endl;

  // Get parameters from configuration
  const eckit::LocalConfiguration saberCentralBlock(conf, "saber central block");
  SaberCentralBlockParametersWrapper saberCentralBlockParamWrapper;
  saberCentralBlockParamWrapper.validateAndDeserialize(saberCentralBlock);
  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamWrapper.saberCentralBlockParameters;

  // Create dummy FieldSet (for xb and fg)
  atlas::FieldSet dummyFs;

  // Create dummy time
  util::DateTime dummyTime(1977, 5, 25, 0, 0, 0);

  // Define central variables
  oops::Variables centralVars = incVars;

  // Get active variables
  oops::Variables activeVars =
    saberCentralBlockParams.activeVars.value().get_value_or(centralVars);

  // Initialize vector of FieldSet
  std::vector<atlas::FieldSet> fsetVec;

  // Read input fields (on model increment geometry)
  std::vector<eckit::LocalConfiguration> inputFieldConfs;
  inputFieldConfs = saberCentralBlockParams.inputFieldConfs.value().get_value_or(inputFieldConfs);
  readInputFields(geom,
                  activeVars,
                  dummyTime,
                  inputFieldConfs,
                  fsetVec);

  // Create central block
  saberCentralBlock_.reset(SaberCentralBlockFactory::create(
                           geom.generic(),
                           geom.variableSizes(activeVars),
                           centralVars,
                           saberCentralBlockParams,
                           dummyFs,
                           dummyFs,
                           fsetVec,
                           geom.timeComm().rank()));

  // Check that active variables are present in central variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(centralVars.has(var));
  }

  oops::Log::trace() << "Localization:Localization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Localization<MODEL>::~Localization() {
  oops::Log::trace() << "Localization:~Localization destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::randomize(Increment_ & dx) const {
  oops::Log::trace() << "Localization:randomize starting" << std::endl;

  // Random vector (necessary for some SABER blocks)
  dx.random();

  // Central block randomization
  saberCentralBlock_->randomize(dx.fieldSet());

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "Localization:randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  oops::Log::trace() << "Localization:multiply starting" << std::endl;

  // Central block multiplication
  saberCentralBlock_->multiply(dx.fieldSet());

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "Localization:multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::print(std::ostream & os) const {
  os << "Localization:print not implemeted yet";
}

// -----------------------------------------------------------------------------

}  // namespace saber
