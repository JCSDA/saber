/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberCentralBlockBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include <boost/noncopyable.hpp>

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/oops/SaberCentralBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberCentralBlockBase::SaberCentralBlockBase(const eckit::Configuration & conf)
  : name_(conf.getString("saber block name")) {}

// -----------------------------------------------------------------------------

SaberCentralBlockFactory::SaberCentralBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberCentralBlockFactory."
                       << std::endl;
    ABORT("Element already registered in saber::SaberCentralBlockFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

SaberCentralBlockBase * SaberCentralBlockFactory::create(
  const eckit::mpi::Comm & comm,
  const oops::GeometryData & geometryData,
  const std::vector<size_t> & activeVariableSizes,
  const eckit::Configuration & conf,
  const atlas::FieldSet & xb,
  const atlas::FieldSet & fg,
  const std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << "SaberCentralBlockBase::create starting" << std::endl;
  const std::string id = conf.getString("saber block name");
  typename std::map<std::string, SaberCentralBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberCentralBlockFactory." << std::endl;
    ABORT("Element does not exist in saber::SaberCentralBlockFactory.");
  }
  SaberCentralBlockBase * ptr = jsb->second->make(comm, geometryData, activeVariableSizes, conf,
                                                  xb, fg, fsetVec);
  oops::Log::trace() << "SaberCentralBlockBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberCentralBlockParametersBase>
SaberCentralBlockFactory::createParameters(const std::string &name) {
  typename std::map<std::string, SaberCentralBlockFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in saber::SaberCentralBlockFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace saber
