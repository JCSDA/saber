/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/covariance/GSI_Covariance.interface.h"
#include "saber/gsi/grid/GSI_Grid.h"
#include "saber/oops/SaberBlockBase.h"


using atlas::option::levels;
using atlas::option::name;

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class CovarianceParameters : public GridParameters {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters, GridParameters)
};

// -------------------------------------------------------------------------------------------------

template <typename MODEL>
class Covariance : public SaberBlockBase<MODEL> {
  typedef oops::Geometry<MODEL> Geometry_;
  typedef oops::State<MODEL>    State_;

 public:
  static const std::string classname() {return "saber::gsi::Covariance";}

  typedef CovarianceParameters Parameters_;

  Covariance(const Geometry_ &, const Parameters_ &, const State_ &, const State_ &);
  virtual ~Covariance();

  void randomize(atlas::FieldSet *) const override;
  void multiply(atlas::FieldSet *) const override;
  void inverseMultiply(atlas::FieldSet *) const override;
  void multiplyAD(atlas::FieldSet *) const override;
  void inverseMultiplyAD(atlas::FieldSet *) const override;

 private:
  void print(std::ostream &) const override;
  // Fortran LinkedList key
  CovarianceKey keySelf_;
  // Variables
  std::vector<std::string> variables_;
  // Function space
  std::unique_ptr<atlas::functionspace::PointCloud> gsiGridFuncSpace_;
  // Grid
  Grid grid_;
};

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
Covariance<MODEL>::Covariance(const Geometry_ & geom, const Parameters_ & params,
                                      const State_ & xbg, const State_ & xfg)
  : SaberBlockBase<MODEL>(params), variables_(), gsiGridFuncSpace_(), grid_(geom.getComm(), params)
{
  oops::Log::trace() << classname() << "::Covariance starting" << std::endl;
  util::Timer timer(classname(), "Covariance");

  // Assert that there is no variable change in this block
  ASSERT(params.inputVars.value() == params.outputVars.value());
  variables_ = params.inputVars.value().variables();

  // Function space
  gsiGridFuncSpace_.reset(new atlas::functionspace::PointCloud(grid_.functionSpace()->get()));

  // Need to convert background and first guess to Atlas and GSI grid.

  // Create covariance module
  gsi_covariance_create_f90(keySelf_, geom.getComm(), params.toConfiguration());

  oops::Log::trace() << classname() << "::Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
Covariance<MODEL>::~Covariance() {
  oops::Log::trace() << classname() << "::~Covariance starting" << std::endl;
  util::Timer timer(classname(), "~Covariance");
  gsi_covariance_delete_f90(keySelf_);
  oops::Log::trace() << classname() << "::~Covariance done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::randomize(atlas::FieldSet * saberFields) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  util::Timer timer(classname(), "randomize");

  // Ignore incoming fields and create new ones based on the block function space
  // ----------------------------------------------------------------------------
  atlas::FieldSet newFields = atlas::FieldSet();

  // Loop over saber (model) fields and create corresponding fields on gsi grid
  for (auto sabField : *saberFields) {
      // Get the name
      const auto fieldName = name(sabField.name());

      // Ensure that the field name is in the input/output list
      const std::string fieldNameStr = fieldName.getString("name");
      if (std::find(variables_.begin(), variables_.end(), fieldNameStr) == variables_.end()) {
        ABORT("Field " + fieldNameStr + " not found in the " + classname() + " variables.");
      }

      // Create the gsi grid field and add to Fieldset
      newFields.add(gsiGridFuncSpace_->createField<double>(fieldName | levels(sabField.levels())));
  }

  // Replace whatever fields are coming in with the gsi grid fields
  *saberFields = newFields;

  // Call implementation
  gsi_covariance_randomize_f90(keySelf_, saberFields->get());
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::multiply(atlas::FieldSet * saberFields) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  util::Timer timer(classname(), "multiply");
  gsi_covariance_multiply_f90(keySelf_, saberFields->get());
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::inverseMultiply(atlas::FieldSet * saberFields) const {
  ABORT(classname() + "inverseMultiply: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::multiplyAD(atlas::FieldSet * saberFields) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  util::Timer timer(classname(), "multiplyAD");
  gsi_covariance_multiply_ad_f90(keySelf_, saberFields->get());
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::inverseMultiplyAD(atlas::FieldSet * saberFields) const {
  ABORT(classname() + "inverseMultiplyAD: not implemented");
}

// -------------------------------------------------------------------------------------------------

template<typename MODEL>
void Covariance<MODEL>::print(std::ostream & os) const {
  os << classname();
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
