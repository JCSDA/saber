/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_
#define SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/generic/instantiateLocalizationFactory.h"
#include "oops/generic/LocalizationGeneric.h"

#include "saber/oops/LocalizationBUMP.h"
#include "saber/oops/LocalizationID.h"

namespace saber {

template <typename MODEL> void instantiateLocalizationFactory() {
  oops::instantiateLocalizationFactory<MODEL>();
  static oops::LocalizationGenericMaker<MODEL, LocalizationBUMP<MODEL> >
    makerBUMP_("BUMP");
  static oops::LocalizationGenericMaker<MODEL, LocalizationID<MODEL> >
    makerID_("ID");
}

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATELOCALIZATIONFACTORY_H_
