/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/blocks/SaberBlockParametersBase.h"

namespace oops {
  class FieldSet3D;
}

namespace saber {

// -----------------------------------------------------------------------------

class SaberOuterBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberOuterBlockBase(const SaberBlockParametersBase & params,
                               const util::DateTime & validTime)
    : blockName_(params.saberBlockName), skipInverse_(params.skipInverse),
      filterMode_(params.filterMode), validTime_(validTime) {}
  virtual ~SaberOuterBlockBase() {}

  // Accessor

  // To inner Geometry data
  virtual const oops::GeometryData & innerGeometryData() const = 0;

  // To inner variables
  virtual const oops::Variables & innerVars() const = 0;

  // Application methods

  // Block multiplication
  virtual void multiply(oops::FieldSet3D &) const = 0;

  // Block multiplication adjoint
  virtual void multiplyAD(oops::FieldSet3D &) const = 0;

  // Block left inverse multiplication
  virtual void leftInverseMultiply(oops::FieldSet3D &) const

    {throw eckit::NotImplemented("leftInverseMultiply not implemented yet for the block "
      + blockName_, Here());}

  // Setup / calibration methods

  // Read block data
  virtual void read()
    {throw eckit::NotImplemented("read not implemented yet for the block " + this->blockName());}

  // Read model fields
  virtual std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const
    {return {};}
  virtual void setReadFields(const std::vector<oops::FieldSet3D> &) {}

  // Direct calibration
  virtual void directCalibration(const oops::FieldSets &)
    {throw eckit::NotImplemented("directCalibration not implemented yet for the block "
      + this->blockName(), Here());}

  // Iterative calibration
  virtual void iterativeCalibrationInit()
    {throw eckit::NotImplemented("iterativeCalibrationInit not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void iterativeCalibrationUpdate(const oops::FieldSet3D &)
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void iterativeCalibrationFinal()
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + this->blockName(), Here());}

  // Dual resolution setup
  virtual void dualResolutionSetup(const oops::GeometryData &)
    {throw eckit::NotImplemented("dualResolutionSetup not implemented yet for the block "
      + this->blockName(), Here());}

  // Write block data
  virtual void write() const {}

  // Write model fields
  virtual std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
     {return {};}

  // Generate inner FieldSet (for the inverse test)
  virtual oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
                                                 const oops::Variables & innerVars) const
    {return oops::randomFieldSet3D(validTime_,
                                   innerGeometryData.comm(),
                                   innerGeometryData.functionSpace(),
                                   innerVars);}

  // Generate outer FieldSet (for the inverse test)
  virtual oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
                                                 const oops::Variables & outerVars) const
    {return oops::randomFieldSet3D(validTime_,
                                   outerGeometryData.comm(),
                                   outerGeometryData.functionSpace(),
                                   outerVars);}

  // Compare FieldSets (for the inverse test)
  virtual bool compareFieldSets(const oops::FieldSet3D & fset3D1,
                                const oops::FieldSet3D & fset3D2,
                                const double & tol) const
    {return fset3D1.compare_with(fset3D2, tol, util::ToleranceType::normalized_absolute);}

  // Non-virtual methods

  // Return block name
  const std::string blockName() const {return blockName_;}

  // Return flag to skip inverse application
  const bool skipInverse() const {return skipInverse_;}

  // Return flag to skip inverse application
  const bool filterMode() const {return filterMode_;}

  // Return date/time
  const util::DateTime validTime() const {return validTime_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::Variables &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &,
             const oops::Variables &) const;

  // Adjoint test
  void adjointTest(const oops::GeometryData &,
                   const oops::Variables &,
                   const oops::GeometryData &,
                   const oops::Variables &,
                   const double &) const;

  // Left-inverse test
  void inverseTest(const oops::GeometryData &,
                   const oops::Variables &,
                   const oops::GeometryData &,
                   const oops::Variables &,
                   const oops::Variables &,
                   const oops::Variables &,
                   const double &,
                   const double &) const;

 private:
  const std::string blockName_;
  const bool skipInverse_;
  const bool filterMode_;
  const util::DateTime validTime_;
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

class SaberOuterBlockFactory;

// -----------------------------------------------------------------------------

class SaberOuterBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberOuterBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberOuterBlockFactory>
    saberOuterBlockParameters{"saber block name", this};
};

// -----------------------------------------------------------------------------

class SaberOuterBlockFactory {
 public:
  static std::unique_ptr<SaberOuterBlockBase> create(const oops::GeometryData &,
                                                     const oops::Variables &,
                                                     const eckit::Configuration &,
                                                     const SaberBlockParametersBase &,
                                                     const oops::FieldSet3D &,
                                                     const oops::FieldSet3D &);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberOuterBlockFactory() = default;

 protected:
  explicit SaberOuterBlockFactory(const std::string &name);

 private:
  virtual std::unique_ptr<SaberOuterBlockBase> make(const oops::GeometryData &,
                                                    const oops::Variables &,
                                                    const eckit::Configuration &,
                                                    const SaberBlockParametersBase &,
                                                    const oops::FieldSet3D &,
                                                    const oops::FieldSet3D &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberOuterBlockFactory * > & getMakers() {
    static std::map < std::string, SaberOuterBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberOuterBlockMaker : public SaberOuterBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<SaberOuterBlockBase> make(const oops::GeometryData & outerGeometryData,
                                            const oops::Variables & outerVars,
                                            const eckit::Configuration & covarConf,
                                            const SaberBlockParametersBase & params,
                                            const oops::FieldSet3D & xb,
                                            const oops::FieldSet3D & fg) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<T>(outerGeometryData, outerVars,
                               covarConf, stronglyTypedParams, xb, fg);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberOuterBlockMaker(const std::string & name) : SaberOuterBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterBlockBase::read(const oops::Geometry<MODEL> & geom,
                               const oops::Variables & vars) {
  oops::Log::trace() << "SaberOuterBlockBase::read starting" << std::endl;

  // Read fieldsets as increments
  std::vector<oops::FieldSet3D> fsetVec;
  for (const auto & input : this->getReadConfs()) {
    // Create increment
    oops::Increment<MODEL> dx(geom, vars, validTime_);
    dx.read(input.second);
    oops::Log::test() << "Norm of input parameter " << input.first
                      << ": " << dx.norm() << std::endl;
    fsetVec.push_back(dx.fieldSet());
    fsetVec.back().name() = input.first;
  }
  this->setReadFields(fsetVec);

  oops::Log::trace() << "SaberOuterBlockBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberOuterBlockBase::write(const oops::Geometry<MODEL> & geom,
                                const oops::Variables & vars) const {
  oops::Log::trace() << "SaberOuterBlockBase::write starting" << std::endl;

  // Get vector of configuration/FieldSet pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs
    = this->fieldsToWrite();

  // Create increment
  oops::Increment<MODEL> dx(geom, vars, validTime_);

  // Loop and write
  for (const auto & output : outputs) {
    dx.fromFieldSet(output.second.fieldSet());
    oops::Log::test() << "Norm of output parameter " << output.second.name()
                      << ": " << dx.norm() << std::endl;
    dx.write(output.first);
  }

  oops::Log::trace() << "SaberOuterBlockBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
