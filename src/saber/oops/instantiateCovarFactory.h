/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_INSTANTIATECOVARFACTORY_H_
#define SABER_OOPS_INSTANTIATECOVARFACTORY_H_

#include "oops/base/instantiateCovarFactory.h"

#include "saber/oops/ErrorCovariance.h"
#include "saber/oops/instantiateLocalizationFactory.h"
#include "saber/oops/instantiateSaberBlockFactory.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateCovarFactory() {
  oops::instantiateCovarFactory<MODEL>();
  static oops::CovarMaker<MODEL, ErrorCovariance<MODEL> > makerSABER_("SABER");

  instantiateLocalizationFactory<MODEL>();
  instantiateSaberBlockFactory();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATECOVARFACTORY_H_
