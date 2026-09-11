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

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/GeometryData.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/blocks/BlockParametersBase.h"

// Forward declaration
namespace oops {
  template <typename MODEL> class Geometry;
  template <typename MODEL> class Increment;
  class FieldSets;
}

namespace saber {

// -----------------------------------------------------------------------------

class CentralBlockBase : public util::Printable,
                              private eckit::NonCopyable {
 public:
  explicit CentralBlockBase(const BlockParametersBase & params,
                            const util::DateTime & validTime,
                            const oops::GeometryData & geometryData,
                            const oops::Variables & centralVars)
    : validTime_(validTime),
      blockName_(params.saberBlockName),
      geometryData_(geometryData),
      centralVars_(centralVars) {}
  virtual ~CentralBlockBase() {}

  // Application methods

  // Block randomization
  virtual void randomize(oops::FieldSet3D &) const;

  // Block multiplication
  virtual void multiply(oops::FieldSet3D &) const;

  // Return the diagonal variance fieldset of this central block. Default
  // throws so every central block must declare it (correlation-only blocks
  // return a unit-variance fieldset; StdDev-style blocks return sigma^2).
  virtual oops::FieldSet3D variance() const
    {throw eckit::NotImplemented("variance not implemented yet for the block "
      + blockName_, Here());}

  // Setup / calibration methods

  // Read block data
  virtual void read()
    {throw eckit::NotImplemented("read not implemented yet for the block " + blockName_,
      Here());}

  // Read model files
  virtual std::vector<std::pair<std::string, eckit::LocalConfiguration>> getReadConfs() const
    {return {};}
  virtual void setReadFields(const std::vector<oops::FieldSet3D> &) {}


  // Direct calibration
  virtual void directCalibration(const oops::FieldSets &)
    {throw eckit::NotImplemented("directCalibration not implemented yet for the block "
      + blockName_, Here());}

  // Iterative calibration
  virtual void iterativeCalibrationInit()
    {throw eckit::NotImplemented("iterativeCalibrationInit not implemented yet for the block "
      + blockName_, Here());}
  virtual void iterativeCalibrationUpdate(const oops::FieldSet3D &)
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + blockName_, Here());}
  virtual void iterativeCalibrationFinal()
    {throw eckit::NotImplemented("iterativeCalibrationUpdate not implemented yet for the block "
      + blockName_, Here());}

  // Write block data
  virtual void write() const {}

  // Write model files
  virtual std::vector<std::pair<eckit::LocalConfiguration, oops::FieldSet3D>> fieldsToWrite() const
     {return {};}

  // Square-root formulation
  virtual size_t ctlVecSize() const
    {throw eckit::NotImplemented("ctlVecSize not implemented yet for the block "
      + blockName_, Here());}
  virtual void randomCtlVec(atlas::Field &, const size_t &) const;
  virtual void multiplySqrt(const atlas::Field &, oops::FieldSet3D &, const size_t &) const
    {throw eckit::NotImplemented("multiplySqrt not implemented yet for the block "
      + blockName_, Here());}
  virtual void multiplySqrtAD(const oops::FieldSet3D &, atlas::Field &, const size_t &) const
    {throw eckit::NotImplemented("multiplySqrtAD not implemented yet for the block "
      + blockName_, Here());}

  // Non-virtual methods

  // Return block name
  const std::string blockName() const {return blockName_;}

  // Return date/time
  const util::DateTime validTime() const {return validTime_;}

  // Return geometry data
  const oops::GeometryData & geometryData() const {return geometryData_;}

  // Return central variables
  const oops::Variables & centralVars() const {return centralVars_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::Variables &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &) const;

 protected:
  const util::DateTime validTime_;

 private:
  const std::string blockName_;
  const oops::GeometryData & geometryData_;
  const oops::Variables centralVars_;

  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

class CentralBlockFactory;

// -----------------------------------------------------------------------------

class CentralBlockParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CentralBlockParametersWrapper, Parameters)
 public:
  oops::RequiredPolymorphicParameter<BlockParametersBase, CentralBlockFactory>
    saberCentralBlockParameters{"saber block name", this};

  const BlockParametersBase & blockParams() const
    {return this->saberCentralBlockParameters;
  }
};

// -----------------------------------------------------------------------------

class CentralBlockFactory {
 public:
  static std::unique_ptr<CentralBlockBase> create(const oops::GeometryData &,
                                                  const oops::Variables &,
                                                  const eckit::Configuration &,
                                                  const BlockParametersBase &,
                                                  const oops::FieldSet3D &,
                                                  const oops::FieldSet3D &);

  static std::unique_ptr<BlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~CentralBlockFactory() = default;

 protected:
  explicit CentralBlockFactory(const std::string &name);

 private:
  virtual std::unique_ptr<CentralBlockBase> make(const oops::GeometryData &,
                                                 const oops::Variables &,
                                                 const eckit::Configuration &,
                                                 const BlockParametersBase &,
                                                 const oops::FieldSet3D &,
                                                 const oops::FieldSet3D &) = 0;

  virtual std::unique_ptr<BlockParametersBase> makeParameters() const = 0;

  static std::map < std::string, CentralBlockFactory * > & getMakers() {
    static std::map < std::string, CentralBlockFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class CentralBlockMaker : public CentralBlockFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<CentralBlockBase> make(const oops::GeometryData & geometryData,
                                         const oops::Variables & outerVars,
                                         const eckit::Configuration & covarConf,
                                         const BlockParametersBase & params,
                                         const oops::FieldSet3D & xb,
                                         const oops::FieldSet3D & fg) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::make_unique<T>(geometryData, outerVars, covarConf,
                               stronglyTypedParams, xb, fg);
  }

  std::unique_ptr<BlockParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit CentralBlockMaker(const std::string & name) : CentralBlockFactory(name) {}
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void CentralBlockBase::read(const oops::Geometry<MODEL> & geom,
                                 const oops::Variables & vars) {
  oops::Log::trace() << "CentralBlockBase::read starting" << std::endl;

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

  oops::Log::trace() << "CentralBlockBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CentralBlockBase::write(const oops::Geometry<MODEL> & geom) const {
  oops::Log::trace() << "CentralBlockBase::write starting" << std::endl;

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

  oops::Log::trace() << "CentralBlockBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
