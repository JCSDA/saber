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

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/gsi/covariance/GSI_Covariance.interface.h"
#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

class CovarianceParameters : public SaberCentralBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(CovarianceParameters, SaberCentralBlockParametersBase)

 public:
  // File containing grid and coefficients
  oops::RequiredParameter<std::string> GSIFile{"gsi error covariance file", this};
  oops::RequiredParameter<std::string> GSINML{"gsi berror namelist file", this};
  oops::RequiredParameter<std::string> GSIVGRD{"gsi akbk", this};

  // Handle vertical top-2-bottom and vice-verse wrt to GSI
  oops::Parameter<bool> vflip{"flip vertical grid", true, this};

  // Processor layout
  oops::Parameter<size_t> layoutx{"processor layout x direction", 1, this};
  oops::Parameter<size_t> layouty{"processor layout y direction", 1, this};

  // Debugging mode
  oops::Parameter<bool> debugMode{"debugging mode", false, this};
  oops::Parameter<bool> bypassGSI{"debugging bypass gsi", false, this};
  oops::Parameter<bool> bypassGSIbe{"debugging deep bypass gsi B error", false, this};
};

// -------------------------------------------------------------------------------------------------

class Covariance : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::gsi::Covariance";}

  typedef CovarianceParameters Parameters_;

  Covariance(const eckit::mpi::Comm &,
     const atlas::FunctionSpace &,
     const atlas::FieldSet &,
     const std::vector<size_t> &,
     const eckit::Configuration &,
     const atlas::FieldSet &,
     const atlas::FieldSet &,
     const std::vector<atlas::FieldSet> &);
  virtual ~Covariance();

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

 private:
  void print(std::ostream &) const override;
  // Fortran LinkedList key
  CovarianceKey keySelf_;
  // Variables
  std::vector<std::string> variables_;
  // GSI grid FunctionSpace
  atlas::FunctionSpace gsiGridFuncSpace_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
