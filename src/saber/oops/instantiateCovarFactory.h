/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SABER_OOPS_INSTANTIATECOVARFACTORY_H_
#define SABER_OOPS_INSTANTIATECOVARFACTORY_H_

#include "oops/base/instantiateCovarFactory.h"

#include "saber/oops/ErrorCovarianceBUMP.h"
#include "saber/oops/ErrorCovarianceGSIRF.h"
#include "saber/oops/ErrorCovarianceID.h"
#include "saber/oops/instantiateLocalizationFactory.h"

namespace saber {

template <typename MODEL> void instantiateCovarFactory() {
  oops::instantiateCovarFactory<MODEL>();
  static oops::CovarMaker<MODEL, ErrorCovarianceBUMP<MODEL> > makerBUMP_("BUMP");
  static oops::CovarMaker<MODEL, ErrorCovarianceGSIRF<MODEL> > makerGSIRF_("GSIRF");
  static oops::CovarMaker<MODEL, ErrorCovarianceID<MODEL> > makerID_("ID");

  saber::instantiateLocalizationFactory<MODEL>();
}

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATECOVARFACTORY_H_
