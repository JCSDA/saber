/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_
#define SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/generic/LocalizationBase.h"

#include "saber/oops/instantiateSaberBlockFactory.h"
#include "saber/oops/Localization.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateLocalizationFactory() {
  static oops::LocalizationMaker<MODEL, Localization<MODEL> >  makerSABER_("SABER");

  instantiateSaberBlockFactory<MODEL>();
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_
