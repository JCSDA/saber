/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/interpolation/GSI_Interpolation.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberOuterBlockBase.h"


namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberOuterBlockMaker<gsi::Interpolation> 
  makerGSI_Interpolation_("gsi interpolation to model grid");

// -------------------------------------------------------------------------------------------------

Interpolation::Interpolation(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf), outputFunctionSpace_(outputFunctionSpace)
{
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;
  util::Timer timer(classname(), "Interpolation");

  // Deserialize configuration
  InterpolationParameters params;
  params.deserialize(conf);

  // Setup GSI grid
  grid_ = Grid(comm, params.grid.value());

  // Input geometry and variables
  inputFunctionSpace_ = atlas::functionspace::PointCloud(grid_.functionSpace());
  inputExtraFields_ = outputExtraFields; // TODO: interpolate that?
  inputVars_ = params.outputVars.value();

  // Object wide copy of the variables
  variables_ = inputVars_.variables();

  // Get number of levels
  gsiLevels_ = grid_.levels();

  // Create the interpolator
  interpolator_.reset(new UnstructuredInterpolation(conf,
                                                    inputFunctionSpace_,
                                                    outputFunctionSpace_,
                                                    nullptr, comm));

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

Interpolation::~Interpolation() {
  oops::Log::trace() << classname() << "::~Interpolation starting" << std::endl;
  util::Timer timer(classname(), "~Interpolation");
  oops::Log::trace() << classname() << "::~Interpolation done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");

  // Create empty Model fieldset
  atlas::FieldSet modFields = atlas::FieldSet();

  // Loop over saber (gsi) fields and create corresponding model fields
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = atlas::option::name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the model field and add to Fieldset
      modFields.add(outputFunctionSpace_.createField<double>(atlas::option::name(fieldName)
      | atlas::option::levels(sabField.levels())));
  }

  // Do the interpolation from GSI grid to model grid
  interpolator_->apply(fset, modFields);

  // Replace the saber (model) fields with the GSI fields
  fset = modFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");

  // Create empty GSI fieldset
  atlas::FieldSet gsiFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding GSI fields
  for (auto sabField : fset) {
      // Get the name
      const auto fieldName = atlas::option::name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the field and add to Fieldset
      gsiFields.add(inputFunctionSpace_.createField<double>(atlas::option::name(fieldName)
        | atlas::option::levels(sabField.levels())));
  }

  // Do the adjoint of interpolation from GSI grid to model grid
  interpolator_->apply_ad(fset, gsiFields);

  // Replace the saber (model) fields with the GSI fields
  fset = gsiFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void Interpolation::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::info() << classname()
                    << "::calibrationInverseMultiply not meaningful so fieldset unchanged"
                    << std::endl;

}

// -------------------------------------------------------------------------------------------------

void Interpolation::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
