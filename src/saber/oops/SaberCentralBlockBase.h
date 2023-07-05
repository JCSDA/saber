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

// TODO(AS): only need a forward declaration of Geometry (and Increment?)
#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "saber/oops/SaberBlockParametersBase.h"

// Forward declaration
namespace saber {
  class SaberBlockChain;
}

namespace saber {

// -----------------------------------------------------------------------------

class SaberCentralBlockBase : public util::Printable, private boost::noncopyable {
 public:
  explicit SaberCentralBlockBase(const SaberBlockParametersBase & params)
    : blockName_(params.saberBlockName) {}
  virtual ~SaberCentralBlockBase() {}

  // Application methods

  // Block multiplication
  virtual void randomize(atlas::FieldSet &) const = 0;

  // Block randomization
  virtual void multiply(atlas::FieldSet &) const = 0;

  // Setup / calibration methods

  // Set localization blockchain
  virtual void setLocalization(std::unique_ptr<SaberBlockChain>)
    {ABORT("setLocalization not implemented yet for the block " + this->blockName());}

  // Read block data
  virtual void read()
    {ABORT("read not implemented yet for the block " + this->blockName());}

  // Read model files
  virtual std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToRead()
    {return {};}

  // Apply inflation on ensemble members
  virtual void applyInflation(std::vector<atlas::FieldSet> &)
    {ABORT("applyInflation not implemented yet for the block " + this->blockName());}

  // Direct calibration
  virtual void directCalibration(const std::vector<atlas::FieldSet> &)
    {ABORT("directCalibration not implemented yet for the block " + this->blockName());}

  // Iterative calibration
  virtual void iterativeCalibrationInit()
    {ABORT("iterativeCalibrationInit not implemented yet for the block " + this->blockName());}
  virtual void iterativeCalibrationUpdate(const atlas::FieldSet &)
    {ABORT("iterativeCalibrationUpdate not implemented yet for the block " + this->blockName());}
  virtual void iterativeCalibrationFinal()
    {ABORT("iterativeCalibrationUpdate not implemented yet for the block " + this->blockName());}

  // Dual resolution setup
  virtual void dualResolutionSetup(const oops::GeometryData &)
    {ABORT("dualResolutionSetup not implemented yet for the block " + this->blockName());}

  // Write block data
  virtual void write() const {}

  // Write model files
  virtual std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> fieldsToWrite() const
     {return {};}

  // Non-virtual methods

  // Return block name
  std::string blockName() const {return blockName_;}

  // Read model fields
  template <typename MODEL>
  void read(const oops::Geometry<MODEL> &,
            const oops::Variables &,
            const util::DateTime &);

  // Write model fields
  template <typename MODEL>
  void write(const oops::Geometry<MODEL> &,
             const oops::Variables &,
             const util::DateTime &) const;

  // Adjoint test
  void adjointTest(const eckit::mpi::Comm &,
                   const oops::GeometryData &,
                   const oops::Variables &,
                   const double & adjointTolerance = 1.0e-12) const;

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
  static SaberCentralBlockBase * create(const oops::GeometryData &,
                                        const oops::Variables &,
                                        const eckit::Configuration &,
                                        const SaberBlockParametersBase &,
                                        const atlas::FieldSet &,
                                        const atlas::FieldSet &,
                                        const util::DateTime &,
                                        const size_t & timeRank = 0);

  static std::unique_ptr<SaberBlockParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SaberCentralBlockFactory() = default;

 protected:
  explicit SaberCentralBlockFactory(const std::string &name);

 private:
  virtual SaberCentralBlockBase * make(const oops::GeometryData &,
                                       const oops::Variables &,
                                       const eckit::Configuration &,
                                       const SaberBlockParametersBase &,
                                       const atlas::FieldSet &,
                                       const atlas::FieldSet &,
                                       const util::DateTime &,
                                       const size_t &) = 0;

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

  SaberCentralBlockBase * make(const oops::GeometryData & geometryData,
                               const oops::Variables & outerVars,
                               const eckit::Configuration & covarConf,
                               const SaberBlockParametersBase & params,
                               const atlas::FieldSet & xb,
                               const atlas::FieldSet & fg,
                               const util::DateTime & validTime,
                               const size_t & timeRank) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(geometryData, outerVars, covarConf,
                 stronglyTypedParams, xb, fg, validTime, timeRank);
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
                                 const oops::Variables & vars,
                                 const util::DateTime & date) {
  oops::Log::trace() << "SaberCentralBlockBase::read starting" << std::endl;

  // Read fieldsets as increments
  for (auto & input : this->fieldsToRead()) {
    // Create increment
    oops::Increment<MODEL> dx(geom, vars, date);
    dx.read(input.first);
    oops::Log::test() << "Norm of input parameter " << input.second.name()
                      << ": " << dx.norm() << std::endl;
    util::copyFieldSet(dx.fieldSet(), input.second);
  }

  oops::Log::trace() << "SaberCentralBlockBase::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralBlockBase::write(const oops::Geometry<MODEL> & geom,
                                  const oops::Variables & vars,
                                  const util::DateTime & date) const {
  oops::Log::trace() << "SaberCentralBlockBase::write starting" << std::endl;

  // Get vector of FieldSet/configuration pairs
  std::vector<std::pair<eckit::LocalConfiguration, atlas::FieldSet>> outputs
    = this->fieldsToWrite();

  // Create increment
  oops::Increment<MODEL> dx(geom, vars, date);

  // Write fieldsets as increments
  for (const auto & output : outputs) {
    dx.fieldSet() = util::copyFieldSet(output.second);
    dx.synchronizeFields();
    oops::Log::test() << "Norm of output parameter " << output.second.name()
                      << ": " << dx.norm() << std::endl;
    dx.write(output.first);
  }

  oops::Log::trace() << "SaberCentralBlockBase::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
