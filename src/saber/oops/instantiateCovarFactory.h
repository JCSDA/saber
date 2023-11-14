/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/instantiateCovarFactory.h"

#include "saber/blocks/instantiateBlockChainFactory.h"
#include "saber/oops/ErrorCovariance.h"
#include "saber/oops/instantiateLocalizationFactory.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateCovarFactory() {
  oops::instantiateCovarFactory<MODEL>();
  static oops::CovarMaker<MODEL, ErrorCovariance<MODEL> > makerSABER_("SABER");
  instantiateLocalizationFactory<MODEL>();
  instantiateBlockChainFactory<MODEL>();
}

// -----------------------------------------------------------------------------

}  // namespace saber
