/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/covariance/Covariance.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/gsi/covariance/Covariance.interface.h"
#include "saber/gsi/grid/Grid.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberCentralBlockMaker<StaticCovariance> makerStaticCovariance_("gsi static covariance");
static SaberCentralBlockMaker<HybridCovariance> makerHybridCovariance_("gsi hybrid covariance");

// -------------------------------------------------------------------------------------------------

StaticCovariance::StaticCovariance(const oops::GeometryData & geometryData,
                       const oops::Variables & centralVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()),
    params_(params), variables_(params.activeVars.value().get_value_or(centralVars)),
    gsiGridFuncSpace_(geometryData.functionSpace()), comm_(&geometryData.comm()),
    xb_(xb.validTime(), xb.commGeom()), fg_(fg.validTime(), fg.commGeom()),
    validTime_(xb.validTime())
{
  oops::Log::trace() << classname() << "::StaticCovariance starting" << std::endl;
  xb_.shallowCopy(xb);
  fg_.shallowCopy(fg);
  oops::Log::trace() << classname() << "::StaticCovariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

StaticCovariance::~StaticCovariance() {
  oops::Log::trace() << classname() << "::~StaticCovariance starting" << std::endl;
  util::Timer timer(classname(), "~StaticCovariance");
  gsi_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~StaticCovariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void StaticCovariance::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on gsi grid
  for (const auto & sabField : fset) {
      // Ensure that the field name is in the variables list
      if (!variables_.has(sabField.name())) {
        throw eckit::Exception("Field " + sabField.name() + " not found in the " + classname() +
          " variables.", Here());
      }

      // Create the gsi grid field and add to Fieldset
      newFields.add(gsiGridFuncSpace_.createField<double>(atlas::option::name(sabField.name())
        | atlas::option::levels(sabField.shape(1))));
  }

  // Replace whatever fields are coming in with the gsi grid fields
  fset.fieldSet() = newFields;

  // Call implementation
  gsi_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void StaticCovariance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  std::vector<const atlas::field::FieldSetImpl*> fsetptrs(1);
  fsetptrs[0] = fset.get();
  gsi_covariance_multiply_f90(keySelf_, fsetptrs.size(), fsetptrs.data());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void StaticCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;
  // Create covariance module
  std::vector<const atlas::field::FieldSetImpl*> fsetXbptrs(1);
  std::vector<const atlas::field::FieldSetImpl*> fsetFgptrs(1);
  std::vector<const util::DateTime *> timesptrs(1);
  fsetXbptrs[0] = xb_.get();
  fsetFgptrs[0] = fg_.get();
  timesptrs[0]  = &validTime_;
  gsi_covariance_create_f90(keySelf_, *comm_, params_.readParams.value()->toConfiguration(),
                            fsetXbptrs.size(), fsetXbptrs.data(), fsetFgptrs.data(),
                            timesptrs.data());
  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void StaticCovariance::print(std::ostream & os) const {
  os << classname();
}
// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
