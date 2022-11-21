/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/State.h"
#include "oops/util/FieldSetOperations.h"

#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralBlockBase.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberCentralTBlockParameters : public SaberCentralBlockParametersWrapper {
  OOPS_CONCRETE_PARAMETERS(SaberCentralTBlockParameters, SaberCentralBlockParametersWrapper)

 public:
  typedef oops::StateEnsembleParameters<MODEL>               StateEnsembleParameters_;
  typedef oops::IncrementEnsembleFromStatesParameters<MODEL> IncrementEnsembleFromStatesParameters_;
  typedef oops::IncrementEnsembleParameters<MODEL>           IncrementEnsembleParameters_;

  /// Ensemble parameters
  oops::OptionalParameter<IncrementEnsembleFromStatesParameters_> ensemble{"ensemble", this};

  /// Ensemble perturbations parameters
  oops::OptionalParameter<IncrementEnsembleParameters_> ensemblePert{"ensemble pert", this};

  /// Ensemble base parameters
  oops::OptionalParameter<StateEnsembleParameters_> ensembleBase{"ensemble base", this};
  /// Ensemble state parameters for the ensemble pairs that would be subtracted from the base
  /// ensemble
  oops::OptionalParameter<StateEnsembleParameters_> ensemblePairs{"ensemble pairs", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class SaberCentralTBlock {
  typedef oops::Geometry<MODEL>                              Geometry_;
  typedef oops::Increment<MODEL>                             Increment_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::StateEnsembleParameters<MODEL>               StateEnsembleParameters_;
  typedef oops::IncrementEnsembleFromStatesParameters<MODEL> IncrementEnsembleFromStatesParameters_;
  typedef oops::IncrementEnsembleParameters<MODEL>           IncrementEnsembleParameters_;
  typedef oops::IncrementEnsemble<MODEL>                     Ensemble_;

 public:
  typedef SaberCentralTBlockParameters<MODEL> Parameters_;

  static const std::string classname() {return "saber::SaberCentralTBlock<MODEL>";}

  SaberCentralTBlock(const Geometry_ &,
                     const oops::GeometryData &,
                     const oops::Variables &,
                     const Parameters_ &,
                     const State_ &,
                     const State_ &,
                     const double & adjointTolerance = -1.0);
  ~SaberCentralTBlock() {}

  void randomize(atlas::FieldSet &) const;
  void multiply(atlas::FieldSet &) const;

 private:
  void print(std::ostream &);
  std::unique_ptr<SaberCentralBlockBase> saberCentralBlock_;
  std::unique_ptr<Ensemble_> ensemble_;
};

// ----------------------------------------------------------------------------

template <typename MODEL>
SaberCentralTBlock<MODEL>::SaberCentralTBlock(const Geometry_ & geom,
                                              const oops::GeometryData & geometryData,
                                              const oops::Variables & outerVars,
                                              const Parameters_ & params,
                                              const State_ & xb,
                                              const State_ & fg,
                                              const double & adjointTolerance)
{
  oops::Log::trace() << classname() << "::SaberCentralTBlock starting" << std::endl;
  // Check for 4D matrices (not ready yet)
  if (geom.timeComm().size() > 1) {
    ABORT("SaberCentralTBlock not ready for 4D matrices yet");
  }

  // Create central block
  const SaberBlockParametersBase & saberCentralBlockParams =
    params.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  // Get active variables
  oops::Variables activeVars = saberCentralBlockParams.activeVars.value().get_value_or(outerVars);

  // Check that active variables are present in outer variables
  for (const auto & var : activeVars.variables()) {
    ASSERT(outerVars.has(var));
  }

  // Read input fields (on model increment geometry)
  std::vector<eckit::LocalConfiguration> inputFields;
  inputFields = saberCentralBlockParams.inputFields.value().get_value_or(inputFields);
  std::vector<atlas::FieldSet> fsetVec = readInputFields(geom,
                                                         activeVars,
                                                         xb.validTime(),
                                                         inputFields);

  // Create central block
  saberCentralBlock_.reset(SaberCentralBlockFactory::create(
                           geometryData,
                           geom.variableSizes(activeVars),
                           activeVars,
                           saberCentralBlockParams,
                           xb.fieldSet(),
                           fg.fieldSet(),
                           fsetVec));

  // Load ensemble
  // TODO(Benjamin): assumption here: ensemble geometry is the increment geometry
  const boost::optional<IncrementEnsembleFromStatesParameters_>
    &ensemble = params.ensemble.value();
  const boost::optional<IncrementEnsembleParameters_> &ensemblePert = params.ensemblePert.value();
  const boost::optional<StateEnsembleParameters_> &ensembleBase = params.ensembleBase.value();
  const boost::optional<StateEnsembleParameters_> &ensemblePairs = params.ensemblePairs.value();

  if (ensemble != boost::none) {
    // Ensemble of states, perturbation using the mean
    oops::Log::info() << "Info     : Ensemble of states, perturbation using the mean"
                      << std::endl;
    State_ xx(geom, activeVars, xb.validTime());
    ensemble_.reset(new Ensemble_(*ensemble, xx, xx, geom, activeVars));
  } else if (ensemblePert) {
    // Increment ensemble from increments on disk
    oops::Log::info() << "Info     : Increment ensemble from increments on disk" << std::endl;
    ensemble_.reset(new Ensemble_(geom, activeVars, *ensemblePert));
  } else if ((ensembleBase != boost::none) &&
             (ensemblePairs != boost::none)) {
    // Increment ensemble from difference of two states
     oops::Log::info() << "Info     : Increment ensemble from difference of two states"
                       << std::endl;
     ensemble_.reset(new Ensemble_(geom, activeVars, *ensembleBase, *ensemblePairs));
  }

  if (adjointTolerance >= 0.0) {
    // Adjoint test

    // Variables sizes
    std::vector<size_t> outerVariableSizes = geom.variableSizes(outerVars);

    // Create random inner FieldSet
    atlas::FieldSet innerFset =
      util::createRandomFieldSet(geometryData.functionSpace(),
                           outerVariableSizes,
                           outerVars);

    // Copy inner FieldSet
    atlas::FieldSet innerFsetSave = util::copyFieldSet(innerFset);

    // Create random outer FieldSet
    atlas::FieldSet outerFset =
      util::createRandomFieldSet(geometryData.functionSpace(),
                           outerVariableSizes,
                           outerVars);

    // Copy outer FieldSet
    atlas::FieldSet outerFsetSave = util::copyFieldSet(outerFset);

    // Apply forward multiplication only (self-adjointness test)
    saberCentralBlock_->multiply(innerFset);
    saberCentralBlock_->multiply(outerFset);

    // Compute adjoint test
    const double dp1 = util::dotProductFieldSets(innerFset, outerFsetSave, activeVars,
                                                 geom.getComm());
    const double dp2 = util::dotProductFieldSets(outerFset, innerFsetSave, activeVars,
                                                 geom.getComm());

    oops::Log::info() << "Info     : Adjoint test for central block "
                      << saberCentralBlockParams.saberBlockName.value()
                      << ": y^t (Ax) = " << dp1 << ": x^t (Ay) = " << dp2 << std::endl;
    ASSERT(abs(dp1) > 0.0);
    ASSERT(abs(dp2) > 0.0);
    oops::Log::test() << "Adjoint test for central block "
                      << saberCentralBlockParams.saberBlockName.value();
    if (0.5*abs(dp1-dp2)/(dp1+dp2) < adjointTolerance) {
      oops::Log::test() << " passed" << std::endl;
    } else {
      oops::Log::test() << " failed" << std::endl;
      ABORT("Adjoint test failure");
    }
  }

  oops::Log::trace() << classname() << "::SaberCentralTBlock done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralTBlock<MODEL>::randomize(atlas::FieldSet & fset) const
{
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  if (ensemble_) {
    // Initialization
    util::zeroFieldSet(fset);

    // Loop over ensemble members
    for (unsigned int ie = 0; ie < ensemble_->size(); ++ie) {
      // Temporary copy for this ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(fset);

      // Randomize SABER central block
      saberCentralBlock_->randomize(fsetMem);

      // Schur product
      util::multiplyFieldSets(fsetMem, (*ensemble_)[ie].fieldSet());

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }

    // Normalize result
    const double rk = 1.0/sqrt(static_cast<double>(ensemble_->size()-1));
    util::multiplyFieldSet(fset, rk);
  } else {
    // Randomize SABER central block
    saberCentralBlock_->randomize(fset);
  }

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void SaberCentralTBlock<MODEL>::multiply(atlas::FieldSet & fset) const
{
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (ensemble_) {
    // Initialization
    atlas::FieldSet fsetInit = util::copyFieldSet(fset);
    util::zeroFieldSet(fset);

    // Loop over ensemble members
    for (unsigned int ie = 0; ie < ensemble_->size(); ++ie) {
      // Temporary copy for this ensemble member
      atlas::FieldSet fsetMem = util::copyFieldSet(fsetInit);

      // First schur product
      util::multiplyFieldSets(fsetMem, (*ensemble_)[ie].fieldSet());

      // Apply localization
      saberCentralBlock_->multiply(fsetMem);

      // Second schur product
      util::multiplyFieldSets(fsetMem, (*ensemble_)[ie].fieldSet());

      // Add up member contribution
      util::addFieldSets(fset, fsetMem);
    }

    // Normalize result
    const double rk = 1.0/static_cast<double>(ensemble_->size()-1);
    util::multiplyFieldSet(fset, rk);
  } else {
    // Apply SABER central block
    saberCentralBlock_->multiply(fset);
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
