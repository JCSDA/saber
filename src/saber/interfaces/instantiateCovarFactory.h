/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SABER_INTERFACE_INSTANTIATECOVARFACTORY_H_
#define SABER_INTERFACE_INSTANTIATECOVARFACTORY_H_

#include "oops/base/instantiateCovarFactory.h"

#include "saber/interfaces/ErrorCovariance4DBUMP.h"
#include "saber/interfaces/ErrorCovarianceBUMP.h"
#include "saber/interfaces/ErrorCovarianceGSIRF.h"
#include "saber/interfaces/instantiateLocalizationFactory.h"

namespace saber {

template <typename MODEL> void instantiateCovarFactory() {
  oops::instantiateCovarFactory<MODEL>();
  static oops::CovarMaker<MODEL, ErrorCovarianceBUMP<MODEL> > makerBUMP_("BUMP");
  static oops::CovarMaker<MODEL, ErrorCovarianceGSIRF<MODEL> > makerGSIRF_("GSIRF");
  static oops::Covar4DMaker<MODEL, ErrorCovariance4DBUMP<MODEL> >  makerBUMP4D_("BUMP");

  oops::instantiateLocalizationFactory<MODEL>();
}

}  // namespace saber

#endif  // SABER_INTERFACE_INSTANTIATECOVARFACTORY_H_
