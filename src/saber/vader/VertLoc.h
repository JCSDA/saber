/*
 * (C) Crown Copyright 2023 Met Office
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

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberOuterBlockBase.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------
class VertLocCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VertLocCovarianceParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> locFileName{"localization matrix file name", this};
  oops::RequiredParameter<std::string> locFieldName{"localization field name in file", this};
  oops::OptionalParameter<std::string> covStatFileName{"pressure file name", this};
  oops::OptionalParameter<std::string> meanPressFieldName{"pressure field name in pressure file",
                                                          this};
};

class VertLocParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(VertLocParameters, SaberBlockParametersBase)
 public:
  oops::RequiredParameter<VertLocCovarianceParameters>
    VertLocParams{"localization data", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
  oops::RequiredParameter<int> nlevels{"number of vertical levels", this};
  oops::RequiredParameter<int> truncation{"number of vertical modes", this};
};

// -----------------------------------------------------------------------------
/// \brief  This saber block applies a vertical localisation to active variables
class VertLoc : public SaberOuterBlockBase {
 public:
  static const std::string classname() {return "saber::vader::VertLoc";}

  typedef VertLocParameters Parameters_;

  VertLoc(const oops::GeometryData &,
          const oops::Variables &,
          const eckit::Configuration &,
          const Parameters_ &,
          const atlas::FieldSet &,
          const atlas::FieldSet &,
          const util::DateTime &);
  virtual ~VertLoc();

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(atlas::FieldSet &) const override;
  void multiplyAD(atlas::FieldSet &) const override;
  void leftInverseMultiply(atlas::FieldSet &) const override;

 private:
  void readLocMat(const std::string filepath, const std::string fieldname,
                  atlas::array::ArrayView<double, 2>& fview);
  void readPressVec(const std::string filepath, const std::string fieldname,
                           atlas::array::ArrayView<double, 1> & fview);
  void print(std::ostream &) const override;
  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
  oops::Variables activeVars_;
  int nlevs;
  int nmods;
  std::string ncfilepath_;
  std::string locfilepath_;
  std::string locFieldName_;
  std::string meanPressFieldName_;
  Eigen::MatrixXd Umatrix_;
};

// -----------------------------------------------------------------------------

}  // namespace vader
}  // namespace saber
