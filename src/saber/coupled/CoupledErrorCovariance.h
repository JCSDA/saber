/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/coupled/TraitCoupled.h"
#include "oops/coupled/UtilsCoupled.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/coupled/CoupledErrorCovarianceParameters.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------
/// \brief Coupled error covariance.
/// \details Implements block-diagonal error covariance for coupled models, with
/// optional common outer blocks.
/// If common outer blocks are specified, two identical common outer blocks chains
/// are created, one for each component model, operating on its geometry.
/// oops::GlobalInterpolator is used to interpolate backgrounds and inputs to the
/// component models' geometries.
template <typename MODEL1, typename MODEL2>
class CoupledErrorCovariance
  : public oops::ModelSpaceCovarianceBase<oops::TraitCoupled<MODEL1, MODEL2>>,
    public util::Printable,
    private util::ObjectCounter<CoupledErrorCovariance<MODEL1, MODEL2>> {
  typedef oops::TraitCoupled<MODEL1, MODEL2>                   COUPLED;
  typedef oops::Geometry<COUPLED>                              Geometry_;
  typedef oops::Increment<COUPLED>                             Increment_;
  typedef oops::Increment4D<COUPLED>                           Increment4D_;
  typedef oops::State4D<COUPLED>                               State4D_;

 public:
  static const std::string classname() {return "saber::CoupledErrorCovariance";}

  CoupledErrorCovariance(const Geometry_ &, const oops::Variables &,
                  const eckit::Configuration &,
                  const State4D_ &, const State4D_ &);
  virtual ~CoupledErrorCovariance() = default;

  void multiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {this->doMultiply(dxi, dxo);}

 private:
  CoupledErrorCovariance(const CoupledErrorCovariance&);
  CoupledErrorCovariance& operator=(const CoupledErrorCovariance&);

  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  void print(std::ostream &) const override;

  /// Chain of blocks (hybrid or ensemble or parametric) for each component model
  std::unique_ptr<SaberBlockChainBase> blockChain1_;
  std::unique_ptr<SaberBlockChainBase> blockChain2_;
  /// @brief  Interpolation operators to the other model's geometry (if needed)
  std::unique_ptr<oops::GlobalInterpolator> interp12_;
  std::unique_ptr<oops::GlobalInterpolator> interp21_;
  /// Variables for each component model
  const std::vector<oops::Variables> incVars_;
  /// @brief  Common outer block chains, one for each model (if needed)
  std::unique_ptr<SaberOuterBlockChain> commonOuterBlockChain1_;
  std::unique_ptr<SaberOuterBlockChain> commonOuterBlockChain2_;

  /// Helper method to create block chain for a component model
  template<typename MODEL>
  std::unique_ptr<SaberBlockChainBase> createBlockChain(
    const oops::Geometry<MODEL> & geom,
    const oops::Variables & vars,
    const ErrorCovarianceParameters & params,
    const oops::State4D<MODEL> & xb,
    const oops::State4D<MODEL> & fg);
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
CoupledErrorCovariance<MODEL1, MODEL2>::CoupledErrorCovariance(const Geometry_ & geom,
                                                                const oops::Variables & incVars,
                                                                const eckit::Configuration & config,
                                                                const State4D_ & xb,
                                                                const State4D_ & fg)
  : oops::ModelSpaceCovarianceBase<COUPLED>(geom, config, xb, fg),
    incVars_(oops::splitVariables(incVars, geom.geometry().variables()))
{
  oops::Log::trace() << "CoupledErrorCovariance::CoupledErrorCovariance starting" << std::endl;
  util::Timer timer(classname(), "CoupledErrorCovariance");
  // Deserialize parameters and fill configuration with missing values
  CoupledErrorCovarianceParameters<MODEL1, MODEL2> params;
  params.deserialize(config);
  eckit::LocalConfiguration fullConf;
  params.serialize(fullConf);

  // Change xb/fg to inner loop resolution
  State4D_ xb_lowres(geom, xb);
  State4D_ fg_lowres(geom, fg);
  // Create State4D objects for each component model
  oops::State4D<MODEL1> xb1 = share_state1(xb_lowres);
  oops::State4D<MODEL2> xb2 = share_state2(xb_lowres);
  oops::State4D<MODEL1> fg1 = share_state1(fg_lowres);
  oops::State4D<MODEL2> fg2 = share_state2(fg_lowres);

  // Create block chains for each component model
  blockChain1_ = createBlockChain<MODEL1>(geom.geometry().geometry1(), incVars_[0],
                                          params.errorCov1.value(), xb1, fg1);
  blockChain2_ = createBlockChain<MODEL2>(geom.geometry().geometry2(), incVars_[1],
                                          params.errorCov2.value(), xb2, fg2);

  // Create interpolation operators to other component's geometry if needed
  if (params.interp.value() && params.commonOuterBlocks.value()) {
    // Check that there are no overlapping variables between the two components
    oops::Variables overlap = xb1.variables();
    overlap.intersection(xb2.variables());
    ASSERT(overlap.size() == 0);
    overlap = fg1.variables();
    overlap.intersection(fg2.variables());
    ASSERT(overlap.size() == 0);

    // Create interpolators
    interp12_ = std::make_unique<oops::GlobalInterpolator>(
             params.interp.value()->interp1.value(),
             geom.geometry().geometry1().generic(),
             geom.geometry().geometry2().generic().functionSpace(),
             geom.geometry().geometry2().generic().comm());
    interp21_ = std::make_unique<oops::GlobalInterpolator>(
             params.interp.value()->interp2.value(),
             geom.geometry().geometry2().generic(),
             geom.geometry().geometry1().generic().functionSpace(),
             geom.geometry().geometry1().generic().comm());

    // Local copy of background and first guess that can undergo interpolation
    oops::FieldSet4D fset4dXb1(xb.times(), xb.commTime(), geom.geometry().getComm());
    oops::FieldSet4D fset4dFg1(fg.times(), fg.commTime(), geom.geometry().getComm());
    oops::FieldSet4D fset4dXb2(xb.times(), xb.commTime(), geom.geometry().getComm());
    oops::FieldSet4D fset4dFg2(fg.times(), fg.commTime(), geom.geometry().getComm());

    for (size_t jt = 0; jt < xb.size(); ++jt) {
      // Get background and first guess fieldsets for both components
      oops::FieldSet3D fset1xb = oops::copyFieldSet3D(xb1[jt].fieldSet());
      oops::FieldSet3D fset2xb = oops::copyFieldSet3D(xb2[jt].fieldSet());
      oops::FieldSet3D fset1fg = oops::copyFieldSet3D(fg1[jt].fieldSet());
      oops::FieldSet3D fset2fg = oops::copyFieldSet3D(fg2[jt].fieldSet());
      // Interpolate states to the other component geometry
      atlas::FieldSet fset1xb_interp, fset2xb_interp;
      atlas::FieldSet fset1fg_interp, fset2fg_interp;
      interp12_->apply(fset1xb.fieldSet(), fset1xb_interp);
      interp12_->apply(fset1fg.fieldSet(), fset1fg_interp);
      interp21_->apply(fset2xb.fieldSet(), fset2xb_interp);
      interp21_->apply(fset2fg.fieldSet(), fset2fg_interp);
      // Add interpolated fields to the relevant fieldsets
      for (auto & field : fset2xb_interp) fset1xb.add(field);
      for (auto & field : fset1xb_interp) fset2xb.add(field);
      for (auto & field : fset2fg_interp) fset1fg.add(field);
      for (auto & field : fset1fg_interp) fset2fg.add(field);
      // Assign to the 4D fieldsets
      fset4dXb1[jt].shallowCopy(fset1xb);
      fset4dXb2[jt].shallowCopy(fset2xb);
      fset4dFg1[jt].shallowCopy(fset1fg);
      fset4dFg2[jt].shallowCopy(fset2fg);
    }
    // Create common outer block chains for each component model
    commonOuterBlockChain1_ = std::make_unique<SaberOuterBlockChain>(
                       geom.geometry().geometry1().generic(),
                       incVars,
                       fset4dXb1,
                       fset4dFg1,
                       fullConf,
                       *params.commonOuterBlocks.value());
    commonOuterBlockChain2_ = std::make_unique<SaberOuterBlockChain>(
                       geom.geometry().geometry2().generic(),
                       incVars,
                       fset4dXb2,
                       fset4dFg2,
                       fullConf,
                       *params.commonOuterBlocks.value());
  }
  oops::Log::trace() << "CoupledErrorCovariance::CoupledErrorCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
template<typename MODEL>
std::unique_ptr<SaberBlockChainBase>
  CoupledErrorCovariance<MODEL1, MODEL2>::createBlockChain(
    const oops::Geometry<MODEL> & geom,
    const oops::Variables & vars,
    const ErrorCovarianceParameters & params,
    const oops::State4D<MODEL> & xb,
    const oops::State4D<MODEL> & fg) {

  // Local copy of background and first guess that can undergo interpolation
  std::unique_ptr<oops::FieldSet4D> fset4dXb;
  std::unique_ptr<oops::FieldSet4D> fset4dFg;

  // Change resolution if needed
  if (params.changeBackgroundResolution) {
    const oops::State4D<MODEL> xb_lowres(geom, xb);
    const oops::State4D<MODEL> fg_lowres(geom, fg);
    const oops::FieldSet4D fset4dXbTmp(xb_lowres);
    const oops::FieldSet4D fset4dFgTmp(fg_lowres);
    fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
    fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
  } else {
    const oops::FieldSet4D fset4dXbTmp(xb);
    const oops::FieldSet4D fset4dFgTmp(fg);
    fset4dXb = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dXbTmp));
    fset4dFg = std::make_unique<oops::FieldSet4D>(oops::copyFieldSet4D(fset4dFgTmp));
  }

  // Initialize outer variables
  const std::vector<std::size_t> vlevs = geom.variableSizes(vars);
  oops::Variables outerVars(vars);
  for (std::size_t i = 0; i < vlevs.size() ; ++i) {
    outerVars[i].setLevels(vlevs[i]);
  }

  // Create block chain
  return SaberBlockChainFactory<MODEL>::create(
       geom,
       outerVars,
       *fset4dXb,
       *fset4dFg,
       params.toConfiguration());
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void CoupledErrorCovariance<MODEL1, MODEL2>::doRandomize(Increment4D_ & dx) const {
  oops::Log::trace() << "CoupledErrorCovariance::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");

  oops::Increment4D<MODEL1> dx1 = share_increment1(dx);
  oops::Increment4D<MODEL2> dx2 = share_increment2(dx);

  // This extra fieldset is only needed for backward compatibility in tests
  oops::FieldSet4D fset4dSum1(dx1.times(), dx1.commTime(), dx1.geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum1.size(); ++jtime) {
    fset4dSum1[jtime].init(blockChain1_->outerFunctionSpace(),
                           blockChain1_->outerVariables(),
                           0.0);
  }

  oops::FieldSet4D fset4dSum2(dx2.times(), dx2.commTime(), dx2.geometry().getComm());
  for (size_t jtime = 0; jtime < fset4dSum2.size(); ++jtime) {
    fset4dSum2[jtime].init(blockChain2_->outerFunctionSpace(),
                           blockChain2_->outerVariables(),
                           0.0);
  }

  // Create FieldSet4D, run randomize on them
  oops::FieldSet4D fset4d1(dx1.times(), dx1.commTime(), dx1.geometry().getComm());
  oops::FieldSet4D fset4d2(dx2.times(), dx2.commTime(), dx2.geometry().getComm());

  blockChain1_->randomize(fset4d1);
  blockChain2_->randomize(fset4d2);

  // If there are common outer blocks, interpolate, apply, and interpolate back
  if (commonOuterBlockChain1_ && commonOuterBlockChain2_) {
    // Interpolate to common geometry
    for (size_t jt = 0; jt < fset4d1.size(); ++jt) {
      atlas::FieldSet fset1_interp, fset2_interp;
      interp12_->apply(fset4d1[jt].fieldSet(), fset1_interp);
      interp21_->apply(fset4d2[jt].fieldSet(), fset2_interp);
      for (auto & field : fset2_interp) {
        fset4d1[jt].add(field);
      }
      for (auto & field : fset1_interp) {
        fset4d2[jt].add(field);
      }
    }
    // Apply common outer blocks
    commonOuterBlockChain1_->applyOuterBlocks(fset4d1);
    commonOuterBlockChain2_->applyOuterBlocks(fset4d2);

    // Remove the variables that were only needed as inputs for the
    // common outer blocks
    fset4d1.removeFields(incVars_[1]);
    fset4d2.removeFields(incVars_[0]);
  }

  // For backward compatibility in tests
  fset4dSum1 += fset4d1;
  fset4dSum2 += fset4d2;

  // ATLAS fieldset to Increment_
  for (size_t jtime = 0; jtime < dx1.size(); ++jtime) {
    dx1[jtime].fromFieldSet(fset4dSum1[jtime].fieldSet());
  }
  for (size_t jtime = 0; jtime < dx2.size(); ++jtime) {
    dx2[jtime].fromFieldSet(fset4dSum2[jtime].fieldSet());
  }

  oops::Log::trace() << "CoupledErrorCovariance::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void CoupledErrorCovariance<MODEL1, MODEL2>::doMultiply(const Increment4D_ & dxi,
                                                         Increment4D_ & dxo) const {
  oops::Log::trace() << "CoupledErrorCovariance::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");

  // Copy input
  dxo = dxi;

  oops::Increment4D<MODEL1> dxo1 = share_increment1(dxo);
  oops::Increment4D<MODEL2> dxo2 = share_increment2(dxo);

  oops::FieldSet4D fset4d1(dxo1);
  oops::FieldSet4D fset4d2(dxo2);

  // Extra fieldset only needed for backward compatibility in tests
  oops::FieldSet4D fset4dSum1 = oops::copyFieldSet4D(fset4d1);
  oops::FieldSet4D fset4dSum2 = oops::copyFieldSet4D(fset4d2);
  fset4dSum1.zero();
  fset4dSum2.zero();

  // If there are common outer blocks, interpolate inputs and apply them
  if (commonOuterBlockChain1_ && commonOuterBlockChain2_) {
    // Interpolate to common geometry
    for (size_t jt = 0; jt < fset4d1.size(); ++jt) {
      atlas::FieldSet fset1_interp, fset2_interp;
      interp12_->apply(fset4d1[jt].fieldSet(), fset1_interp);
      interp21_->apply(fset4d2[jt].fieldSet(), fset2_interp);
      for (auto & field : fset2_interp) {
        fset4d1[jt].add(field);
      }
      for (auto & field : fset1_interp) {
        fset4d2[jt].add(field);
      }
    }

    // Apply common outer blocks
    commonOuterBlockChain1_->applyOuterBlocksAD(fset4d1);
    commonOuterBlockChain2_->applyOuterBlocksAD(fset4d2);

    // Remove the variables that were only needed as inputs for the
    // common outer blocks
    fset4d1.removeFields(incVars_[1]);
    fset4d2.removeFields(incVars_[0]);
  }

  // Apply SABER block chains
  blockChain1_->multiply(fset4d1);
  blockChain2_->multiply(fset4d2);

  // If there are common outer blocks, interpolate, apply, and interpolate back
  if (commonOuterBlockChain1_ && commonOuterBlockChain2_) {
    // Interpolate to common geometry
    for (size_t jt = 0; jt < fset4d1.size(); ++jt) {
      atlas::FieldSet fset1_interp, fset2_interp;
      interp12_->apply(fset4d1[jt].fieldSet(), fset1_interp);
      interp21_->apply(fset4d2[jt].fieldSet(), fset2_interp);
      for (auto & field : fset2_interp) {
        fset4d1[jt].add(field);
      }
      for (auto & field : fset1_interp) {
        fset4d2[jt].add(field);
      }
    }
    // Apply common outer blocks
    commonOuterBlockChain1_->applyOuterBlocks(fset4d1);
    commonOuterBlockChain2_->applyOuterBlocks(fset4d2);

    // Remove the variables that were only needed as inputs for the
    // common outer blocks
    fset4d1.removeFields(incVars_[1]);
    fset4d2.removeFields(incVars_[0]);
  }

  // For backward compatibility in tests
  fset4dSum1 += fset4d1;
  fset4dSum2 += fset4d2;

  // ATLAS fieldset to Increment
  for (size_t jtime = 0; jtime < dxo1.size(); ++jtime) {
    dxo1[jtime].fromFieldSet(fset4dSum1[jtime].fieldSet());
  }
  for (size_t jtime = 0; jtime < dxo2.size(); ++jtime) {
    dxo2[jtime].fromFieldSet(fset4dSum2[jtime].fieldSet());
  }

  oops::Log::trace() << "CoupledErrorCovariance::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void CoupledErrorCovariance<MODEL1, MODEL2>::doInverseMultiply(const Increment4D_ & dxi,
                                                                Increment4D_ & dxo) const {
  oops::Log::trace() << "CoupledErrorCovariance::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");

  // Iterative inverse - this operates on the full coupled increment
  oops::IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

  oops::Log::trace() << "CoupledErrorCovariance::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void CoupledErrorCovariance<MODEL1, MODEL2>::print(std::ostream & os) const {
  oops::Log::trace() << "CoupledErrorCovariance::print starting" << std::endl;
  os << "SABER CoupledErrorCovariance";
  oops::Log::trace() << "CoupledErrorCovariance::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber
