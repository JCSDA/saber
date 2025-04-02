/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <tuple>

#include "saber/blocks/SaberParametricBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {

// Generic constructor for the SaberParametricBlockChain, not templated on
// MODEL. This constructor is only used for localization matrices, when the
// outer geometry cannot be a MODEL geometry.
SaberParametricBlockChain::SaberParametricBlockChain(
                          const oops::GeometryData & outerGeometryData,
                          const oops::Variables & outerVars,
                          oops::FieldSet4D & fset4dXb,
                          oops::FieldSet4D & fset4dFg,
                          const eckit::LocalConfiguration & covarConf,
                          const eckit::Configuration & conf)
  : outerFunctionSpace_(outerGeometryData.functionSpace()),
    outerVariables_(outerVars),
    crossTimeCov_(covarConf.getString("time covariance") == "multivariate duplicated"),
    timeComm_(fset4dXb.commTime()),
    size4D_(fset4dXb.size()) {
  oops::Log::trace() << "SaberParametricBlockChain generic ctor starting" << std::endl;

  // If needed create generic outer block chain
  if (conf.has("saber outer blocks")) {
    std::vector<SaberOuterBlockParametersWrapper> cmpOuterBlocksParams;
    for (const auto & cmpOuterBlockConf : conf.getSubConfigurations("saber outer blocks")) {
      SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
      cmpOuterBlockParamsWrapper.deserialize(cmpOuterBlockConf);
      cmpOuterBlocksParams.push_back(cmpOuterBlockParamsWrapper);
    }
    outerBlockChain_ = std::make_unique<SaberOuterBlockChain>(outerGeometryData,
                                                              outerVariables_,
                                                              fset4dXb,
                                                              fset4dFg,
                                                              covarConf,
                                                              cmpOuterBlocksParams);
  }

  // Set outer geometry data for central block
  const oops::GeometryData & currentOuterGeom = outerBlockChain_ ?
                             outerBlockChain_->innerGeometryData() : outerGeometryData;

  SaberCentralBlockParametersWrapper saberCentralBlockParamsWrapper;
  saberCentralBlockParamsWrapper.deserialize(conf.getSubConfiguration("saber central block"));

  const SaberBlockParametersBase & saberCentralBlockParams =
    saberCentralBlockParamsWrapper.saberCentralBlockParameters;
  oops::Log::info() << "Info     : Creating central block: "
                    << saberCentralBlockParams.saberBlockName.value() << std::endl;

  const auto[currentOuterVars, activeVars]
              = initCentralBlock(currentOuterGeom,
                                 conf,
                                 covarConf,
                                 saberCentralBlockParams,
                                 fset4dXb,
                                 fset4dFg);

  // Check block doesn't expect model fields to be read as this is a generic ctor
  if (centralBlock_->getReadConfs().size() != 0) {
    throw eckit::UserError("The generic constructor of the SABER parametric block chain "
                           "does not allow to read MODEL fields.", Here());
  }

  // Check block doesn't expect calibration, as this could be done with the standard ctor
  if (saberCentralBlockParams.doCalibration()) {
    throw eckit::UserError("The generic constructor of the SABER parametric block chain "
                           "does not allow covariance calibration.", Here());
  }
  if (covarConf.has("dual resolution ensemble configuration")) {
    throw eckit::UserError("The generic constructor of the SABER parametric block chain "
                           "does not allow dual resolution ensemble.", Here());
  }
  if (covarConf.has("output ensemble")) {
    throw eckit::UserError("The generic constructor of the SABER parametric block chain "
                           "does not allow ensemble output.", Here());
  }

  if (saberCentralBlockParams.doRead()) {
    // Read data
    oops::Log::info() << "Info     : Read data" << std::endl;
    centralBlock_->read();

    if (saberCentralBlockParams.forceWrite.value()) {
      // Write data
      oops::Log::info() << "Info     : Write data" << std::endl;
      centralBlock_->write();
    }
  }

  testCentralBlock(covarConf, saberCentralBlockParams, currentOuterGeom, activeVars);

  oops::Log::trace() << "SaberParametricBlockChain generic ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

std::tuple<oops::Variables, oops::Variables>
    SaberParametricBlockChain::initCentralBlock(
        const oops::GeometryData & outerGeom,
        const eckit::Configuration & conf,
        const eckit::LocalConfiguration & covarConf,
        const SaberBlockParametersBase & saberCentralBlockParams,
        const oops::FieldSet4D & fset4dXb,
        const oops::FieldSet4D & fset4dFg) {
  oops::Log::trace() << "SaberParametricBlockChain::initCentralBlock starting" << std::endl;
  // Set outer variables for central block
  const oops::Variables currentOuterVars = outerBlockChain_ ?
                             outerBlockChain_->innerVars() : outerVariables_;

  // Get active variables
  oops::Variables activeVars = getActiveVars(saberCentralBlockParams, currentOuterVars);
  // Check that active variables are present in variables
  for (const auto & var : activeVars) {
    if (!currentOuterVars.has(var)) {
      throw eckit::UserError("Active variable " + var.name() + " is not present in "
                             "outer variables", Here());
    }
  }

  // Create central block
  centralBlock_ = SaberCentralBlockFactory::create(outerGeom,
                                                   activeVars,
                                                   covarConf,
                                                   saberCentralBlockParams,
                                                   fset4dXb[0],
                                                   fset4dFg[0]);

  // Save central function space and variables
  centralFunctionSpace_ = outerGeom.functionSpace();
  centralVars_ = activeVars;

  auto out = std::tuple<oops::Variables, oops::Variables>(currentOuterVars, activeVars);
  oops::Log::trace() << "SaberParametricBlockChain::initCentralBlock exiting..."
                     << std::endl;
  return out;
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::testCentralBlock(
        const eckit::LocalConfiguration & covarConf,
        const SaberBlockParametersBase & saberCentralBlockParams,
        const oops::GeometryData & outerGeom,
        const oops::Variables & activeVars) const {
  oops::Log::trace() << "SaberParametricBlockChain::testCentralBlock starting" << std::endl;
  // Adjoint test
  if (covarConf.getBool("adjoint test")) {
    // Get tolerance
    const double localAdjointTolerance =
      saberCentralBlockParams.adjointTolerance.value().get_value_or(
      covarConf.getDouble("adjoint tolerance"));

    // Run test
    centralBlock_->adjointTest(outerGeom,
                               activeVars,
                               localAdjointTolerance);
  }

  // Square-root test
  if (covarConf.getBool("square-root test")) {
    // Get tolerance
    const double localSqrtTolerance =
      saberCentralBlockParams.sqrtTolerance.value().get_value_or(
      covarConf.getDouble("square-root tolerance"));

    // Run test
    centralBlock_->sqrtTest(outerGeom,
                            activeVars,
                            localSqrtTolerance);
  }
  oops::Log::trace() << "SaberParametricBlockChain::testCentralBlock done" << std::endl;
}

// -----------------------------------------------------------------------------
void SaberParametricBlockChain::filter(oops::FieldSet4D & fset4d) const {
  // Outer blocks for adjoint multiplication or left inverse (acting as filter)
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksFilter(fset4d);
  }

  // No cross-time covariances: apply central block to each of the
  // time slots.
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    centralBlock_->multiply(fset4d[jtime]);
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiply(oops::FieldSet4D & fset4d) const {
  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4d);
  }

  // Central block multiplication
  if (crossTimeCov_) {
    // Duplicated cross-time covariances.
    // Use Mark Buehner's trick to save CPU when applying the same 3D cov/loc for all
    // 3D blocks of the 4D cov/loc matrix:
    // C_4D = ( C_3D C_3D C_3D ) = ( Id ) C_3D ( Id Id Id )
    //        ( C_3D C_3D C_3D )   ( Id )
    //        ( C_3D C_3D C_3D )   ( Id )
    // so if :
    // x_4D = ( x_1 )
    //        ( x_2 )
    //        ( x_3 )
    // then:
    // C_4D x_4D = (Id) C_3D (Id Id Id) (x_1) = (Id) C_3D (x_1+x_2+x_3) = (C_3D ( x_1 + x_2 + x_3 ))
    //             (Id)                 (x_2)   (Id)                      (C_3D ( x_1 + x_2 + x_3 ))
    //             (Id)                 (x_3)   (Id)                      (C_3D ( x_1 + x_2 + x_3 ))
    // Reference in section 3.4.2. of https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2325.
    // Local sum of x1, x2, ...
    for (size_t jtime = 1; jtime < fset4d.size(); ++jtime) {
      fset4d[0] += fset4d[jtime];
    }
    if (timeComm_.rank() > 0) {
      oops::mpi::send(timeComm_, fset4d[0], 0, 0);
    } else {
      // On rank 0 receive other local sums and compute global sum of x1, x2, ...
      oops::FieldSet3D fset3d_tmp = oops::initFieldSet3D(fset4d[0]);
      for (size_t jj = 1; jj < timeComm_.size(); ++jj) {
        oops::mpi::receive(timeComm_, fset3d_tmp, jj, 0);
        fset4d[0] += fset3d_tmp;
      }
      // Compute C * (x1+x2+...)
      centralBlock_->multiply(fset4d[0]);
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].deepCopy(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances: apply central block to each of the
    // time slots.
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->multiply(fset4d[jtime]);
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}


// -----------------------------------------------------------------------------

void SaberParametricBlockChain::randomize(oops::FieldSet4D & fset4d) const {
  // Create central FieldSet4D
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    fset4d[jtime].init(centralFunctionSpace_, centralVars_);
  }

  // Central block randomization
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      centralBlock_->randomize(fset4d[0]);
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].deepCopy(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->randomize(fset4d[jtime]);
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

size_t SaberParametricBlockChain::ctlVecSize() const {
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      // Central block square-root for rank 0
      return centralBlock_->ctlVecSize();
    } else {
      // No control vector
      return 0;
    }
  } else {
    // No cross-time covariances
    return centralBlock_->ctlVecSize()*size4D_;
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiplySqrt(const atlas::Field & cv,
                                             oops::FieldSet4D & fset4d,
                                             const size_t & offset) const {
  // Create central FieldSet4D
  for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
    fset4d[jtime].init(centralFunctionSpace_, centralVars_);
  }

  // Central block square-root
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    if (timeComm_.rank() == 0) {
      // Central block square-root for rank 0
      centralBlock_->multiplySqrt(cv, fset4d[0], offset);
    }
    // Broadcast the result to all tasks
    oops::mpi::broadcast(timeComm_, fset4d[0], 0);
    // Deep copy of the result to all the local time slots
    for (size_t jt = 1; jt < fset4d.local_time_size(); ++jt) {
      fset4d[jt].deepCopy(fset4d[0].fieldSet());
    }
  } else {
    // No cross-time covariances
    size_t index = offset;
    for (size_t jtime = 0; jtime < fset4d.size(); ++jtime) {
      centralBlock_->multiplySqrt(cv, fset4d[jtime], index);
      index += centralBlock_->ctlVecSize();
    }
  }

  // Outer blocks forward multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocks(fset4d);
  }
}

// -----------------------------------------------------------------------------

void SaberParametricBlockChain::multiplySqrtAD(const oops::FieldSet4D & fset4d,
                                               atlas::Field & cv,
                                               const size_t & offset) const {
  // Initialization
  oops::FieldSet4D fset4dCopy = oops::copyFieldSet4D(fset4d);

  // Outer blocks adjoint multiplication
  if (outerBlockChain_) {
    outerBlockChain_->applyOuterBlocksAD(fset4dCopy);
  }

  // Central block square-root adjoint
  if (crossTimeCov_) {
    // Duplicated cross-time covariances
    for (size_t jtime = 1; jtime < fset4dCopy.size(); ++jtime) {
      fset4dCopy[0] += fset4dCopy[jtime];
    }
    if (timeComm_.rank() > 0) {
      oops::mpi::send(timeComm_, fset4dCopy[0], 0, 0);
    } else {
      // On rank 0 receive other local sums and compute global sum of x1, x2, ...
      oops::FieldSet3D fset3d_tmp = oops::initFieldSet3D(fset4dCopy[0]);
      for (size_t jj = 1; jj < timeComm_.size(); ++jj) {
        oops::mpi::receive(timeComm_, fset3d_tmp, jj, 0);
        fset4dCopy[0] += fset3d_tmp;
      }

      // Central block square-root adjoint for rank 0
      centralBlock_->multiplySqrtAD(fset4dCopy[0], cv, offset);
    }
  } else {
    // No cross-time covariances
    size_t index = offset;
    for (size_t jtime = 0; jtime < fset4dCopy.size(); ++jtime) {
      centralBlock_->multiplySqrtAD(fset4dCopy[jtime], cv, index);
      index += centralBlock_->ctlVecSize();
    }
  }
}

// -----------------------------------------------------------------------------


}  // namespace saber
