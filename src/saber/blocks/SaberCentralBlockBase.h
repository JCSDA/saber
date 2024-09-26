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

#include "atlas/field.h"

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/blocks/SaberBlockParametersBase.h"

// Forward declaration
namespace oops {
  template <typename MODEL> class Geometry;
  template <typename MODEL> class Increment;
  class FieldSets;
}

namespace saber {

// -----------------------------------------------------------------------------

class SaberCentralBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberCentralBlockBase(const SaberBlockParametersBase & params,
                                 const util::DateTime & validTime)
    : validTime_(validTime), blockName_(params.saberBlockName) {}
  virtual ~SaberCentralBlockBase() {}

  // Application methods

  // Block multiplication
  virtual void randomize(oops::FieldSet3D &) const = 0;

  // Block randomization
  virtual void multiply(oops::FieldSet3D &) const = 0;

  // Setup / calibration methods

  // Read block data
  virtual void read()
    {throw eckit::NotImplemented("read not implemented yet for the block " + this->blockName(),
      Here());}

  // Read model files
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

  // Write model files
  virtual std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
     {return {};}

  // Square-root formulation
  virtual size_t ctlVecSize() const
    {throw eckit::NotImplemented("ctlVecSize not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const
    {throw eckit::NotImplemented("multiplySqrt not implemented yet for the block "
      + this->blockName(), Here());}
  virtual void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const
    {throw eckit::NotImplemented("multiplySqrtAD not implemented yet for the block "
      + this->blockName(), Here());}

  // Non-virtual methods

  // Return block name
  std::string blockName() const {return blockName_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::Variables &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &) const;

  // Adjoint test
  void adjointTest(const oops::GeometryData &,
                   const oops::Variables &,
                   const double &) const;

  // Square-root test
  void sqrtTest(const oops::GeometryData &,
                const oops::Variables &,
                const double &) const;

 protected:
  const util::DateTime validTime_;

 private:
  std::string blockName_;
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

class SaberCentralBlockFactory;

// -----------------------------------------------------------------------------

class SaberCentralBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SaberCentralBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<SaberBlockParametersBase, SaberCentralBlockFactory>
    saberCentralBlockParameters{"saber block name", this};
};

// -----------------------------------------------------------------------------

class SaberCentralBlockFactory {
 public:
  static std::unique_ptr<SaberCentralBlockBase> create(const oops::GeometryData &,
                                                       const oops::Variables &,
                                                       const eckit::Configuration &,
                                                       const SaberBlockParametersBase &,
                                                       const oops::FieldSet3D &,
                                                       const oops::FieldSet3D &);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberCentralBlockFactory() = default;

 protected:
  explicit SaberCentralBlockFactory(const std::string &name);

 private:
  virtual std::unique_ptr<SaberCentralBlockBase> make(const oops::GeometryData &,
                                                      const oops::Variables &,
                                                      const eckit::Configuration &,
                                                      const SaberBlockParametersBase &,
                                                      const oops::FieldSet3D &,
                                                      const oops::FieldSet3D &) = 0;

  virtual std::unique_ptr<SaberBlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, SaberCentralBlockFactory * > & getMakers() {
    static std::map < std::string, SaberCentralBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class SaberCentralBlockMaker : public SaberCentralBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<SaberCentralBlockBase> make(const oops::GeometryData & geometryData,
                                              const oops::Variables & outerVars,
                                              const eckit::Configuration & covarConf,
                                              const SaberBlockParametersBase & params,
                                              const oops::FieldSet3D & xb,
                                              const oops::FieldSet3D & fg) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<T>(geometryData, outerVars, covarConf,
                               stronglyTypedParams, xb, fg);
  }

  std::unique_ptr<SaberBlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SaberCentralBlockMaker(const std::string & name) : SaberCentralBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralBlockBase::read(const oops::Geometry<MODEL> & geom,
                                 const oops::Variables & vars) {
  oops::Log::trace() << "SaberCentralBlockBase::read starting" << std::endl;

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

  oops::Log::trace() << "SaberCentralBlockBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralBlockBase::write(const oops::Geometry<MODEL> & geom) const {
  oops::Log::trace() << "SaberCentralBlockBase::write starting" << std::endl;

  // Get vector of FieldSet/configuration pairs
  std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> outputs
    = this->fieldsToWrite();

  // Write fieldsets as increments
  for (const auto & output : outputs) {
    oops::Increment<MODEL> dx(geom, output.second.variables(), validTime_);
    dx.fromFieldSet(output.second.fieldSet());
    oops::Log::test() << "Norm of output parameter " << output.second.name()
                      << ": " << dx.norm() << std::endl;
    dx.write(output.first);
  }

  oops::Log::trace() << "SaberCentralBlockBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
