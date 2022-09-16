/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/covariance/GSI_Covariance.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/covariance/GSI_Covariance.interface.h"
#include "saber/gsi/grid/GSI_Grid.h"
#include "saber/gsi/interpolation/GSI_InterpolationImpl.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberCentralBlockMaker<gsi::Covariance> makerGSI_Covariance_("gsi covariance");

// -------------------------------------------------------------------------------------------------

Covariance::Covariance(const eckit::mpi::Comm & comm,
       const atlas::FunctionSpace & functionSpace,
       const atlas::FieldSet & extraFields,
       const std::vector<size_t> & activeVariableSizes,
       const eckit::Configuration & conf,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg,
       const std::vector<atlas::FieldSet> & fsetVec)
  : SaberCentralBlockBase(conf), gsiGridFuncSpace_(functionSpace)
{
  oops::Log::trace() << classname() << "::Covariance starting" << std::endl;
  util::Timer timer(classname(), "Covariance");

  // Deserialize configuration
  CovarianceParameters params;
  params.validateAndDeserialize(conf);

  // Object wide copy of the variables
  variables_ = params.inoutVars.value().variables();

  // Create covariance module
  gsi_covariance_create_f90(keySelf_, comm, params.toConfiguration(), xb.get(), fg.get());

  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Covariance::~Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  gsi_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on gsi grid
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = atlas::option::name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the gsi grid field and add to Fieldset
      newFields.add(gsiGridFuncSpace_.createField<double>(fieldName
        | atlas::option::levels(sabField.levels())));
  }

  // Replace whatever fields are coming in with the gsi grid fields
  fset = newFields;

  // Call implementation
  gsi_covariance_randomize_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  gsi_covariance_multiply_f90(keySelf_, fset.get());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Covariance::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
