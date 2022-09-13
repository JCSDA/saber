/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/oops/SaberBlockBase.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include <boost/noncopyable.hpp>

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

#include "saber/oops/SaberBlockParametersBase.h"

namespace saber {

// -----------------------------------------------------------------------------

SaberBlockBase::SaberBlockBase(const SaberBlockParametersBase & params)
  : iterativeInverse_(params.iterativeInverse.value()), name_(params.saberBlockName.value()) {}

// -----------------------------------------------------------------------------

SaberBlockFactory::SaberBlockFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in saber::SaberBlockFactory." << std::endl;
    ABORT("Element already registered in saber::SaberBlockFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

SaberBlockBase * SaberBlockFactory::create(const atlas::FunctionSpace & functionSpace,
                                           const atlas::FieldSet & extraFields,
                                           const std::vector<size_t> & variableSizes,
                                           const SaberBlockParametersBase & params,
                                           const atlas::FieldSet & xb,
                                           const atlas::FieldSet & fg,
                                           const std::vector<atlas::FieldSet> & fsetVec) {
  oops::Log::trace() << "SaberBlockBase::create starting" << std::endl;
  const std::string id = params.saberBlockName.value();
  typename std::map<std::string, SaberBlockFactory*>::iterator jsb = getMakers().find(id);
  if (jsb == getMakers().end()) {
    oops::Log::error() << id << " does not exist in saber::SaberBlockFactory." << std::endl;
    ABORT("Element does not exist in saber::SaberBlockFactory.");
  }
  SaberBlockBase * ptr = jsb->second->make(functionSpace, extraFields, variableSizes,
                                           params, xb, fg, fsetVec);
  oops::Log::trace() << "SaberBlockBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<SaberBlockParametersBase>
SaberBlockFactory::createParameters(const std::string &name) {
  typename std::map<std::string, SaberBlockFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in saber::SaberBlockFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace saber
