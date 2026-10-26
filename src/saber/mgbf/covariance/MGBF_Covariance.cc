/*
 * (C) Copyright 2024 DOC/NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/mgbf/covariance/MGBF_Covariance.h"

#include <iostream>

#include "atlas/field.h"

#include "saber/mgbf/covariance/MGBF_Covariance.interface.h"
#include "saber/oops/Utilities.h"

namespace saber {
namespace mgbf {

// -------------------------------------------------------------------------------------------------

static CentralBlockMaker<MGBF_Covariance> makerCovariance_("MGBF_covariance");

// -------------------------------------------------------------------------------------------------

MGBF_Covariance::MGBF_Covariance(const oops::GeometryData & geometryData,
                                 const oops::Variables & centralVars,
                                 const eckit::Configuration & covarConf,
                                 const Parameters_ & params,
                                 const oops::FieldSet3D & xb,
                                 const oops::FieldSet3D & fg)
  : CentralBlockBase(params, xb.validTime(), geometryData, centralVars),
     params_(params), variables_(params.activeVars.value().get_value_or(centralVars).variables()),
     mgbfGridFuncSpace_(geometryData.functionSpace()), comm_(&geometryData.comm())
{
  oops::Log::trace() << classname() << "MGBF::Covariance starting" << std::endl;

  util::Timer timer(classname(), "Covariance");
  eckit::LocalConfiguration mgbf_config = params.toConfiguration();
  if (params.doCalibration()) {
    throw eckit::UserError("doCalibration=true is not implemented ", Here());
  }

  // Create covariance module
  mgbf_covariance_create_f90(keySelf_, *comm_, mgbf_config,
                            mgbfGridFuncSpace_.get(), xb.get(), fg.get());

  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

MGBF_Covariance::~MGBF_Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  mgbf_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void MGBF_Covariance::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  for (auto sabField : fset) {
    // Get the name
    const auto fieldName = name(sabField.name());

    // Ensure that the field name is in the input/output list
    const std::string fieldNameStr = fieldName.getString("name");
    if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
      ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
    }
  }

  mgbf_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void MGBF_Covariance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Get member index
  int index_member;
    index_member = 9999;

  // Multiply
  mgbf_covariance_multiply_f90(keySelf_, fset.get(), index_member);

  // Mark all fields as having dirty halos after modification
  for (const auto & fieldname : fset.field_names()) {
    atlas::Field field = fset[fieldname];
    field.set_dirty();  // Mark field as having dirty halos that need to be synchronized
  }

  // Perform the actual halo exchange
  fset.fieldSet().haloExchange();

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void MGBF_Covariance::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace mgbf
}  // namespace saber
