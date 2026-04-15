/*
 * (C) Copyright 2025- UCAR
 * (C) Copyright 2025- NOAA/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/instantiateCovarFactory.h"
#include "oops/coupled/BlockDiagonalCovarianceCoupled.h"
#include "oops/coupled/TraitCoupled.h"

#include "saber/blocks/instantiateBlockChainFactory.h"
#include "saber/coupled/CoupledErrorCovariance.h"
#include "saber/oops/ErrorCovariance.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL1, typename MODEL2> void instantiateCoupledCovarFactory() {
  // covariances both for coupled and single models
  saber::instantiateCovarFactory<oops::TraitCoupled<MODEL1, MODEL2>>();
  saber::instantiateCovarFactory<MODEL1>();
  saber::instantiateCovarFactory<MODEL2>();
  // oops block diagonal covariance
  static oops::CovarMaker<oops::TraitCoupled<MODEL1, MODEL2>,
                          oops::BlockDiagonalCovarianceCoupled<MODEL1, MODEL2> >
                          makerCoupled_("Coupled Block Diagonal");
  // saber coupled covariance
  static oops::CovarMaker<oops::TraitCoupled<MODEL1, MODEL2>,
                          CoupledErrorCovariance<MODEL1, MODEL2> >
                          makerSABERCoupled_("SABER coupled");
}

// -----------------------------------------------------------------------------

}  // namespace saber
