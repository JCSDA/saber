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

#include "saber/oops/SaberCentralTBlock.h"

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
  std::vector<std::reference_wrapper<const oops::GeometryData>> geometryData_;
  std::unique_ptr<SaberCentralTBlock<MODEL>> saberCentralTBlock_;
};

// =============================================================================

template<typename MODEL>
Localization<MODEL>::Localization(const Geometry_ & geom,
                                  const oops::Variables & incVars,
                                  const eckit::Configuration & conf)
  : saberCentralTBlock_()
{
  oops::Log::trace() << "Localization::Localization starting" << std::endl;

  size_t myslot = geom.timeComm().rank();
  if (myslot == 0) {
    // Create dummy time
    util::DateTime dummyTime(1977, 5, 25, 0, 0, 0);

    // Dummy state
    State_ xx(geom, incVars, dummyTime);

    // Get parameters from configuration
    const eckit::LocalConfiguration saberCentralTBlockConf(conf, "saber central block");
    SaberCentralTBlockParameters<MODEL> saberCentralTBlockParams;
    saberCentralTBlockParams.validateAndDeserialize(saberCentralTBlockConf);

    // Create central block wrapper
    saberCentralTBlock_.reset(new SaberCentralTBlock<MODEL>(geom,
                              geom.generic(),
                              incVars,
                              saberCentralTBlockParams,
                              xx,
                              xx));
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
  saberCentralTBlock_->randomize(dx.fieldSet());

  // ATLAS fieldset to Increment_
  dx.synchronizeFields();

  oops::Log::trace() << "Localization:randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Localization<MODEL>::multiply(Increment_ & dx) const {
  oops::Log::trace() << "Localization:multiply starting" << std::endl;

  // Central block multiplication
  saberCentralTBlock_->multiply(dx.fieldSet());

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
